/* *
 * Copyright (c) 2014, James S. Plank and Kevin Greenan
 * All rights reserved.
 *
 * Jerasure - A C/C++ Library for a Variety of Reed-Solomon and RAID-6 Erasure
 * Coding Techniques
 *
 * Revision 2.0: Galois Field backend now links to GF-Complete
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 *  - Neither the name of the University of Tennessee nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
 * WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

/* Jerasure's authors:

   Revision 2.x - 2014: James S. Plank and Kevin M. Greenan
   Revision 1.2 - 2008: James S. Plank, Scott Simmerman and Catherine D. Schuman.
   Revision 1.0 - 2007: James S. Plank
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#include "galois.h"
#include "jerasure.h"
#include <infiniband/verbs_exp.h>
#include <sys/syscall.h>
#include <pthread.h>

#define talloc(type, num) (type *) malloc(sizeof(type)*(num))

static double jerasure_total_xor_bytes = 0;
static double jerasure_total_gf_bytes = 0;
static double jerasure_total_memcpy_bytes = 0;

void jerasure_print_matrix(int *m, int rows, int cols, int w)
{
  int i, j;
  int fw;
  char s[30];
  unsigned int w2;

  if (w == 32) {
    fw = 10;
  } else {
    w2 = (1 << w);
    sprintf(s, "%u", w2-1);
    fw = strlen(s);
  }

  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      if (j != 0) printf(" ");
      printf("%*u", fw, m[i*cols+j]); 
    }
    printf("\n");
  }
}

void jerasure_print_bitmatrix(int *m, int rows, int cols, int w)
{
  int i, j;

  for (i = 0; i < rows; i++) {
    if (i != 0 && i%w == 0) printf("\n");
    for (j = 0; j < cols; j++) {
      if (j != 0 && j%w == 0) printf(" ");
      printf("%d", m[i*cols+j]); 
    }
    printf("\n");
  }
}

int jerasure_make_decoding_matrix(int k, int m, int w, int *matrix, int *erased, int *decoding_matrix, int *dm_ids)
{
  int i, j, *tmpmat;

  j = 0;
  for (i = 0; j < k; i++) {
    if (erased[i] == 0) {
      dm_ids[j] = i;
      j++;
    }
  }

  tmpmat = talloc(int, k*k);
  if (tmpmat == NULL) { return -1; }
  for (i = 0; i < k; i++) {
    if (dm_ids[i] < k) {
      for (j = 0; j < k; j++) tmpmat[i*k+j] = 0;
      tmpmat[i*k+dm_ids[i]] = 1;
    } else {
      for (j = 0; j < k; j++) {
        tmpmat[i*k+j] = matrix[(dm_ids[i]-k)*k+j];
      }
    }
  }

  i = jerasure_invert_matrix(tmpmat, decoding_matrix, k, w);
  free(tmpmat);
  return i;
}

/* Internal Routine */
int jerasure_make_decoding_bitmatrix(int k, int m, int w, int *matrix, int *erased, int *decoding_matrix, int *dm_ids)
{
  int i, j, *tmpmat;
  int index, mindex;

  j = 0;
  for (i = 0; j < k; i++) {
    if (erased[i] == 0) {
      dm_ids[j] = i;
      j++;
    }
  }

  tmpmat = talloc(int, k*k*w*w);
  if (tmpmat == NULL) { return -1; }
  for (i = 0; i < k; i++) {
    if (dm_ids[i] < k) {
      index = i*k*w*w;
      for (j = 0; j < k*w*w; j++) tmpmat[index+j] = 0;
      index = i*k*w*w+dm_ids[i]*w;
      for (j = 0; j < w; j++) {
        tmpmat[index] = 1;
        index += (k*w+1);
      }
    } else {
      index = i*k*w*w;
      mindex = (dm_ids[i]-k)*k*w*w;
      for (j = 0; j < k*w*w; j++) {
        tmpmat[index+j] = matrix[mindex+j];
      }
    }
  }

  i = jerasure_invert_bitmatrix(tmpmat, decoding_matrix, k*w);
  free(tmpmat);
  return i;
}

int jerasure_matrix_decode(int k, int m, int w, int *matrix, int row_k_ones, int *erasures,
                          char **data_ptrs, char **coding_ptrs, int size)
{
  int i, edd, lastdrive;
  int *tmpids;
  int *erased, *decoding_matrix, *dm_ids;

  if (w != 8 && w != 16 && w != 32) return -1;

  erased = jerasure_erasures_to_erased(k, m, erasures);
  if (erased == NULL) return -1;

  /* Find the number of data drives failed */

  lastdrive = k;

  edd = 0;
  for (i = 0; i < k; i++) {
    if (erased[i]) {
      edd++;
      lastdrive = i;
    }
  }
    
  /* You only need to create the decoding matrix in the following cases:

      1. edd > 0 and row_k_ones is false.
      2. edd > 0 and row_k_ones is true and coding device 0 has been erased.
      3. edd > 1

      We're going to use lastdrive to denote when to stop decoding data.
      At this point in the code, it is equal to the last erased data device.
      However, if we can't use the parity row to decode it (i.e. row_k_ones=0
         or erased[k] = 1, we're going to set it to k so that the decoding 
         pass will decode all data.
   */

  if (!row_k_ones || erased[k]) lastdrive = k;

  dm_ids = NULL;
  decoding_matrix = NULL;

  if (edd > 1 || (edd > 0 && (!row_k_ones || erased[k]))) {
    dm_ids = talloc(int, k);
    if (dm_ids == NULL) {
      free(erased);
      return -1;
    }

    decoding_matrix = talloc(int, k*k);
    if (decoding_matrix == NULL) {
      free(erased);
      free(dm_ids);
      return -1;
    }

    if (jerasure_make_decoding_matrix(k, m, w, matrix, erased, decoding_matrix, dm_ids) < 0) {
      free(erased);
      free(dm_ids);
      free(decoding_matrix);
      return -1;
    }
  }

  /* Decode the data drives.  
     If row_k_ones is true and coding device 0 is intact, then only decode edd-1 drives.
     This is done by stopping at lastdrive.
     We test whether edd > 0 so that we can exit the loop early if we're done.
   */

  for (i = 0; edd > 0 && i < lastdrive; i++) {
    if (erased[i]) {
      jerasure_matrix_dotprod(k, w, decoding_matrix+(i*k), dm_ids, i, data_ptrs, coding_ptrs, size);
      edd--;
    }
  }

  /* Then if necessary, decode drive lastdrive */

  if (edd > 0) {
    tmpids = talloc(int, k);
    if (!tmpids) {
      free(erased);
      free(dm_ids);
      free(decoding_matrix);
      return -1;
    }
    for (i = 0; i < k; i++) {
      tmpids[i] = (i < lastdrive) ? i : i+1;
    }
    jerasure_matrix_dotprod(k, w, matrix, tmpids, lastdrive, data_ptrs, coding_ptrs, size);
    free(tmpids);
  }
  
  /* Finally, re-encode any erased coding devices */

  for (i = 0; i < m; i++) {
    if (erased[k+i]) {
      jerasure_matrix_dotprod(k, w, matrix+(i*k), NULL, i+k, data_ptrs, coding_ptrs, size);
    }
  }

  free(erased);
  if (dm_ids != NULL) free(dm_ids);
  if (decoding_matrix != NULL) free(decoding_matrix);

  return 0;
}


int *jerasure_matrix_to_bitmatrix(int k, int m, int w, int *matrix) 
{
  int *bitmatrix;
  int rowelts, rowindex, colindex, elt, i, j, l, x;

  if (matrix == NULL) { return NULL; }

  bitmatrix = talloc(int, k*m*w*w);
  if (!bitmatrix) return NULL;

  rowelts = k * w;
  rowindex = 0;

  for (i = 0; i < m; i++) {
    colindex = rowindex;
    for (j = 0; j < k; j++) {
      elt = matrix[i*k+j];
      for (x = 0; x < w; x++) {
        for (l = 0; l < w; l++) {
          bitmatrix[colindex+x+l*rowelts] = ((elt & (1 << l)) ? 1 : 0);
        }
        elt = galois_single_multiply(elt, 2, w);
      }
      colindex += w;
    }
    rowindex += rowelts * w;
  }
  return bitmatrix;
}

#define err_log     dprintf
#define info_log    dprintf
int g_fd;

struct ibv_device *
find_device(const char *devname)
{
    struct ibv_device **dev_list = NULL;
    struct ibv_device *device = NULL;

    dev_list = ibv_get_device_list(NULL);
    if (!dev_list) {
        err_log(g_fd, "Failed to get IB devices list.\n");
        return NULL;
    }

    if (!devname) {
        device = dev_list[0];
        if (!device)
            err_log(g_fd, "No IB devices found\n");
    } else {
        int i;

        for (i = 0; dev_list[i]; ++i)
            if (!strcmp(ibv_get_device_name(dev_list[i]),
                    devname))
                break;

        device = dev_list[i];
        if (!device)
            err_log(g_fd, "IB device %s not found\n", devname);
    }

/*    ibv_free_device_list(dev_list);*/

    return device;
}

struct ec_context {
    struct ibv_context               *context;
    struct ibv_pd                    *pd;
    struct ibv_exp_ec_calc           *calc;
    struct ibv_exp_ec_calc_init_attr attr;
    int                              block_size;
    struct ibv_sge                   *data_sge;
    struct ibv_sge                   *code_sge;
    struct ibv_exp_ec_mem            mem;
};

typedef enum InitEC_status
{
  NOT_INITIALIZED = 0,
  INIT_SUCCESS,
  INIT_FAILED
} InitEC_status;

#define CONTEXT_NUM 32

struct encoder_context {
  struct ibv_context  *context;
  struct ibv_pd       *pd;
  struct ec_context*  ec_ctx[CONTEXT_NUM];
  unsigned int        cnt;
  /* one global key */
  struct ibv_mr       *mr;
  /* w that device supports */
  uint32_t            ec_w_mask;
  struct ibv_device   *device;
};

struct inargs {
    char    *devname;
    int     k;
    int     m;
    int     w;
    int     *matrix;
    int     size;
};

struct encoder_context    *g_ctx = NULL;
static InitEC_status       g_offload_init_status = NOT_INITIALIZED;

static int
ec_get_sg(struct ibv_sge **data,
          struct ibv_pd *pd,
          int sge_size,
          int block_size,
          char **data_ptrs)
{
    int i;

    *data = calloc(sge_size, sizeof(**data));
    if (!*data) {
        err_log(g_fd, "Failed to allocate data sges\n");
        goto fail;
    }

    for (i = 0; i < sge_size; i++) {
        (*data)[i].lkey = g_ctx->mr->lkey;
        (*data)[i].addr = (uintptr_t)*(data_ptrs + i);
        (*data)[i].length = block_size;
    }

    return 0;

fail:
    return -ENOMEM;
}

static int
alloc_ec_mrs(struct ec_context *ctx, char **data_ptrs, char **code_ptrs)
{
  int err;

  err = ec_get_sg(&ctx->data_sge, ctx->pd, ctx->attr.k, ctx->block_size, data_ptrs);
  if (err)
  return err;

  err = ec_get_sg(&ctx->code_sge, ctx->pd, ctx->attr.m, ctx->block_size, code_ptrs);
  if (err)
  goto free_dbuf;


  ctx->mem.data_blocks = ctx->data_sge;
  ctx->mem.num_data_sge = ctx->attr.k;
  ctx->mem.code_blocks = ctx->code_sge;
  ctx->mem.num_code_sge = ctx->attr.m;
  ctx->mem.block_size = ctx->block_size;
  
  return 0;

free_dbuf:
  free(ctx->data_sge);

  return err;
}

static void free_ec_mr(struct ibv_sge **data)
{
  if (!data || !*data)
    return;

  free(*data);
  *data = NULL;
}

static void
free_ec_mrs(struct ec_context *ctx)
{

    free_ec_mr(&ctx->code_sge);
    free_ec_mr(&ctx->data_sge);
}

void free_ec_ctx(struct ec_context *ctx)
{
  ibv_exp_dealloc_ec_calc(ctx->calc);
  free_ec_mrs(ctx);
  /*free(ctx->encode_matrix);*/
  free(ctx->attr.encode_matrix);
  free(ctx);
}

static int alloc_encode_matrix(int k, int m, int w, int *rs_mat, uint8_t **en_mat)
{
    uint8_t *matrix;
    int i, j;

    matrix = calloc(1, m * k);
    if (!matrix) {
        err_log(g_fd, "Failed to allocate encode matrix\n");
        return -ENOMEM;
    }

    for (i = 0; i < m; i++)
        for (j = 0; j < k; j++)
            matrix[j*m+i] = (uint8_t)rs_mat[i*k+j];

    *en_mat = matrix;

    return 0;
}

struct ec_context *
alloc_ec_ctx2(struct ibv_pd *pd, struct inargs *in)
{
    struct ec_context *ctx;
    int err;

    ctx = calloc(1, sizeof(*ctx));
    if (!ctx) {
        err_log(g_fd, "Failed to allocate EC context\n");
        return NULL;
    }

    ctx->pd = pd;
    ctx->context = pd->context;

    ctx->attr.comp_mask = IBV_EXP_EC_CALC_ATTR_MAX_INFLIGHT |
            IBV_EXP_EC_CALC_ATTR_K |
            IBV_EXP_EC_CALC_ATTR_M |
            IBV_EXP_EC_CALC_ATTR_W |
            IBV_EXP_EC_CALC_ATTR_MAX_DATA_SGE |
            IBV_EXP_EC_CALC_ATTR_MAX_CODE_SGE |
            IBV_EXP_EC_CALC_ATTR_ENCODE_MAT |
            IBV_EXP_EC_CALC_ATTR_AFFINITY |
            IBV_EXP_EC_CALC_ATTR_POLLING;
    ctx->attr.max_inflight_calcs = 1;
    ctx->attr.k = in->k;
    ctx->attr.m = in->m;
    ctx->attr.w = in->w;
    ctx->attr.max_data_sge = in->k;
    ctx->attr.max_code_sge = in->m;
    ctx->attr.affinity_hint = 0;
    ctx->block_size = in->size;

    err = alloc_encode_matrix(ctx->attr.k, ctx->attr.m,
                              ctx->attr.w, in->matrix, &ctx->attr.encode_matrix);
    if (err)
        goto free_mrs;

    ctx->calc = ibv_exp_alloc_ec_calc(ctx->pd, &ctx->attr);

    if (!ctx->calc) {
        err_log(g_fd, "Failed to allocate EC calc\n");
        goto free_mrs;
    }
    
    return ctx;

free_mrs:
    free(ctx);

    return NULL;
}

struct ec_context *
alloc_ec_ctx(struct ibv_pd *pd, struct inargs *in, char **data_ptrs, char **code_ptrs)
{
    struct ec_context *ctx;
    int err;

    ctx = calloc(1, sizeof(*ctx));
    if (!ctx) {
        err_log(g_fd, "Failed to allocate EC context\n");
        return NULL;
    }

    ctx->pd = pd;
    ctx->context = pd->context;

    ctx->attr.comp_mask = IBV_EXP_EC_CALC_ATTR_MAX_INFLIGHT |
            IBV_EXP_EC_CALC_ATTR_K |
            IBV_EXP_EC_CALC_ATTR_M |
            IBV_EXP_EC_CALC_ATTR_W |
            IBV_EXP_EC_CALC_ATTR_MAX_DATA_SGE |
            IBV_EXP_EC_CALC_ATTR_MAX_CODE_SGE |
            IBV_EXP_EC_CALC_ATTR_ENCODE_MAT |
            IBV_EXP_EC_CALC_ATTR_AFFINITY |
            IBV_EXP_EC_CALC_ATTR_POLLING;
    ctx->attr.max_inflight_calcs = 1;
    ctx->attr.k = in->k;
    ctx->attr.m = in->m;
    ctx->attr.w = in->w;
    ctx->attr.max_data_sge = in->k;
    ctx->attr.max_code_sge = in->m;
    ctx->attr.affinity_hint = 0;
    ctx->block_size = in->size;

    err_log(g_fd, "block_size=%d\n", ctx->block_size);

    err = alloc_ec_mrs(ctx, data_ptrs, code_ptrs);
    err_log(g_fd, "after alloc_ec_mrs err %d\n", err);
    if (err)
        goto free_mrs;
    err = alloc_encode_matrix(ctx->attr.k, ctx->attr.m,
                              ctx->attr.w, in->matrix, &ctx->attr.encode_matrix);
    if (err)
        goto free_mrs;

    ctx->calc = ibv_exp_alloc_ec_calc(ctx->pd, &ctx->attr);
    err_log(g_fd, "after ibv_exp_alloc_ec_calc\n");

    if (!ctx->calc) {
        err_log(g_fd, "Failed to allocate EC calc\n");
        goto free_mrs;
    }
    
    return ctx;

free_mrs:
    free_ec_mrs(ctx);
    free(ctx);

    return NULL;
}

void update_ec_ctx(struct ec_context *ctx, char **data_ptrs, char **coding_ptrs)
{
  int i;

  for (i=0; i < ctx->attr.k; i++) {
    ctx->data_sge[i].addr = (uintptr_t)*(data_ptrs + i);
  }
  for (i=0; i < ctx->attr.m; i++) {
    ctx->code_sge[i].addr = (uintptr_t)*(coding_ptrs + i);
  }

}

#define ULLONG_MAX 0xFFFFFFFFFFFFFFFF
void __attribute__ ((constructor)) init()
{
  struct ibv_exp_device_attr dattr;
  int err;

  char fname[100] ={'\0'};
  pid_t pid, tid;

  pid = getpid();
  tid = syscall(SYS_gettid);

  sprintf(fname, "/tmp/debug_ceph.jul_%d_%d", (int)pid, (int)tid);

  g_fd = open(fname, O_RDWR | O_CREAT | O_APPEND , 0666);
  if (g_fd < 0) {
    fprintf(stderr, "ERROR: failed to open file\n");
  }

  dprintf(g_fd, "In init.\n");
  g_offload_init_status = INIT_FAILED;

  struct ibv_device *device;
  device = find_device(NULL);
  if (!device) {
    dprintf(g_fd, "ERROR: init didn't find device\n");
    assert(0);
  }
  err_log(g_fd, "device %s\n", ibv_get_device_name(device));

  g_ctx = calloc(1, sizeof(*g_ctx));
  if (!g_ctx) {
    err_log(g_fd, "Failed to allocate encoder context\n");
    return;
  }

  g_ctx->context = ibv_open_device(device);
  if (!g_ctx->context) {
      err_log(g_fd, "Couldn't get context for %s\n",
                    ibv_get_device_name(device));
      goto free_ctx;
  }
  g_ctx->pd = ibv_alloc_pd(g_ctx->context);
  if (!g_ctx->pd) {
    err_log(g_fd, "Failed to allocate PD\n");
    goto close_device;
  }

  memset(&dattr, 0, sizeof(dattr));
  dattr.comp_mask = IBV_EXP_DEVICE_ATTR_EXP_CAP_FLAGS |
		    IBV_EXP_DEVICE_ATTR_EC_CAPS |
		    IBV_EXP_DEVICE_ATTR_EC_GF_BASE;
  err = ibv_exp_query_device(g_ctx->context, &dattr);
  if (err) {
      err_log(g_fd, "Couldn't query device for EC offload caps. err = %d\n", err);
      goto free_ctx;
  }

  if (!(dattr.exp_device_cap_flags & IBV_EXP_DEVICE_EC_OFFLOAD)) {
      err_log(g_fd, "EC offload not supported by driver.\n");
      goto free_ctx;
  }

  g_ctx->device = device;
  g_ctx->ec_w_mask = dattr.ec_w_mask;

  info_log(g_fd, "EC offload supported by driver.\n");
  info_log(g_fd, "max_ec_calc_inflight_calcs %d\n", dattr.ec_caps.max_ec_calc_inflight_calcs);
  info_log(g_fd, "max_data_vector_count %d\n", dattr.ec_caps.max_ec_data_vector_count);

  struct ibv_exp_reg_mr_in in;
  in.pd = g_ctx->pd;
  in.addr = 0;
  in.length = ULLONG_MAX;
  in.exp_access = IBV_EXP_ACCESS_LOCAL_WRITE | IBV_EXP_ACCESS_ON_DEMAND;
  in.comp_mask = 0;
  in.create_flags = 0;
  g_ctx->mr = ibv_exp_reg_mr(&in);

  if (!g_ctx->mr) {
    err_log(g_fd, "Failed to allocate data MR\n");
    goto close_device;
  }
  /* k=2, m=1, w=8 */
  int *rs_mat; /*reed_sol_vandermonde_coding_matrix(2, 1, 8);*/
  rs_mat = talloc(int, 2);
  if (!rs_mat) {
    err_log(g_fd, "Failed to alloc matrix\n");
    goto close_device;
  }
  rs_mat[0]=rs_mat[1]=1;
  struct inargs alloc_in = {};

  alloc_in.k = 2;
  alloc_in.m = 1;
  alloc_in.w = 8;
  alloc_in.matrix = rs_mat;

  for (int i = 0; i < CONTEXT_NUM; i++) {
    g_ctx->ec_ctx[i] = alloc_ec_ctx2(g_ctx->pd, &alloc_in);
    if (!g_ctx->ec_ctx[i]) {
      err_log(g_fd, "Failed to allocate EC context for i=%d.\n", i);
      goto close_device;
    }
  }
  g_ctx->cnt = 0;
  g_offload_init_status = INIT_SUCCESS;
  return;

close_device:
  ibv_close_device(g_ctx->context);
  /* TBD */
  /* remove allocated contexes that succeeded */
free_ctx:
  free(g_ctx);
  g_ctx = NULL;
}

void __attribute__ ((destructor)) fini()
{
  if (g_offload_init_status != INIT_SUCCESS)
    return;
  int i;
  for (i = 0; i < CONTEXT_NUM; i++)
    free_ec_ctx(g_ctx->ec_ctx[i]);

  ibv_dealloc_pd(g_ctx->pd);

  if (ibv_close_device(g_ctx->context))
      err_log(g_fd, "Couldn't release context\n");

  free(g_ctx);
  err_log(g_fd, "fini: out\n");
  close(g_fd);
}

pthread_mutex_t g_mutex;
void jerasure_matrix_encode(int k, int m, int w, int *matrix,
                          char **data_ptrs, char **coding_ptrs, int size)
{
  int err, i;
  /*pid_t tid;*/
  struct ec_context *choose_ctx = NULL;

  /*tid = syscall(SYS_gettid);*/

  if (w != 8 && w != 16 && w != 32) {
    fprintf(stderr, "ERROR: jerasure_matrix_encode() and w is not 8, 16 or 32\n");
    assert(0);
  }

  if (g_offload_init_status == INIT_FAILED) {
    fprintf(stderr, "In jerasure_matrix_encode, INIT_FAILED\n");
    goto old;
  }

  if (!(g_ctx->ec_w_mask & (1 << (w - 1)))) {
      err_log(g_fd, "W(%d) not supported for given device(%s)\n",
	      w, ibv_get_device_name(g_ctx->context->device));
      goto old;
  }

  struct ibv_exp_ec_mem mem;
  struct ibv_sge data_sge[2] = {};
  struct ibv_sge code_sge[1] = {};

  mem.block_size = size;

  mem.num_data_sge = k;
  /*mem.data_blocks = talloc(struct ibv_sge, k);*/
  mem.data_blocks = data_sge;
  for (i = 0; i < k; i++) {
    /*mem.data_blocks[i].lkey = g_ctx->mr->lkey;
    mem.data_blocks[i].addr = (uintptr_t)*(data_ptrs + i);
    mem.data_blocks[i].length = size;*/
    data_sge[i].lkey = g_ctx->mr->lkey;
    data_sge[i].addr = (uintptr_t)*(data_ptrs + i);
    data_sge[i].length = size;
  }
  

  mem.num_code_sge = m;
  /*mem.code_blocks = talloc(struct ibv_sge, m);*/
  mem.code_blocks = code_sge;
  for (i = 0; i < m; i++) {
   /*mem.code_blocks[i].lkey = g_ctx->mr->lkey;
    mem.code_blocks[i].addr = (uintptr_t)*(coding_ptrs + i);
    mem.code_blocks[i].length = size;*/
    code_sge[i].lkey = g_ctx->mr->lkey;
    code_sge[i].addr = (uintptr_t)*(coding_ptrs + i);
    code_sge[i].length = size;
  }

  unsigned int cnt = 0;

  pthread_mutex_lock(&g_mutex);
  cnt = g_ctx->cnt % CONTEXT_NUM;
  choose_ctx = g_ctx->ec_ctx[cnt];
  g_ctx->cnt++;
  pthread_mutex_unlock(&g_mutex);
  /*err_log(g_fd, "tid %d choose ec_ctx %p (g_ctx->cnt=%u, cnt %u)\n", tid, choose_ctx, g_ctx->cnt, cnt);*/
/* 
  pthread_mutex_lock(&g_mutex);
  if (g_ctx->cnt1 <= g_ctx->cnt2) {
    choose_ctx = g_ctx->ec_ctx;
    err_log(g_fd, "tid %d choose ec_ctx %p cnt1 (cnt1=%u, cnt2 %u)\n", tid, choose_ctx, g_ctx->cnt1, g_ctx->cnt2);
    g_ctx->cnt1 ++;
    choose_cnt = &g_ctx->cnt1;
  } else {
    choose_ctx = g_ctx->ec_ctx2;
    err_log(g_fd, "tid %d choose ec_ctx2 %p cnt2 (cnt1=%u, cnt2 %u)\n", tid, choose_ctx, g_ctx->cnt1, g_ctx->cnt2);
    g_ctx->cnt2 ++;
    choose_cnt = &g_ctx->cnt2;
  }
  pthread_mutex_unlock(&g_mutex);*/

  /*err = ibv_exp_ec_encode_sync(g_ctx->ec_ctx->calc, &(g_ctx->ec_ctx->mem));*/
  err = ibv_exp_ec_encode_sync(choose_ctx->calc, &mem);
  if (err)
    err_log(g_fd, "Failed ibv_exp_ec_encode (%d)\n", err);
 
/*  pthread_mutex_lock(&g_mutex);
  *choose_cnt = *choose_cnt - 1;
  err_log(g_fd, "tid %d unlock choose ctx %p cnt1 (cnt1=%u, cnt2 %u)\n", tid, choose_ctx, g_ctx->cnt1, g_ctx->cnt2);
  pthread_mutex_unlock(&g_mutex);*/
  /*free(mem.data_blocks);
  free(mem.code_blocks);*/
  goto close;

old:
  for (i = 0; i < m; i++) {
    jerasure_matrix_dotprod(k, w, matrix+(i*k), NULL, k+i, data_ptrs, coding_ptrs, size);
  }
close:
  return;
  /*err_log(g_fd, "jerasure_matrix_encode: Out.\n");*/
}

void jerasure_bitmatrix_dotprod(int k, int w, int *bitmatrix_row,
                             int *src_ids, int dest_id,
                             char **data_ptrs, char **coding_ptrs, int size, int packetsize)
{
  int j, sindex, pstarted, index, x, y;
  char *dptr, *pptr, *bdptr, *bpptr;

  if (size%(w*packetsize) != 0) {
    fprintf(stderr, "jerasure_bitmatrix_dotprod - size%c(w*packetsize)) must = 0\n", '%');
    assert(0);
  }

  bpptr = (dest_id < k) ? data_ptrs[dest_id] : coding_ptrs[dest_id-k];

  for (sindex = 0; sindex < size; sindex += (packetsize*w)) {
    index = 0;
    for (j = 0; j < w; j++) {
      pstarted = 0;
      pptr = bpptr + sindex + j*packetsize;
      for (x = 0; x < k; x++) {
        if (src_ids == NULL) {
          bdptr = data_ptrs[x];
        } else if (src_ids[x] < k) {
          bdptr = data_ptrs[src_ids[x]];
        } else {
          bdptr = coding_ptrs[src_ids[x]-k];
        }
        for (y = 0; y < w; y++) {
          if (bitmatrix_row[index]) {
            dptr = bdptr + sindex + y*packetsize;
            if (!pstarted) {
              memcpy(pptr, dptr, packetsize);
              jerasure_total_memcpy_bytes += packetsize;
              pstarted = 1;
            } else {
              galois_region_xor(dptr, pptr, packetsize);
              jerasure_total_xor_bytes += packetsize;
            }
          }
          index++;
        }
      }
    }
  }
}

void jerasure_do_parity(int k, char **data_ptrs, char *parity_ptr, int size) 
{
  int i;

  memcpy(parity_ptr, data_ptrs[0], size);
  jerasure_total_memcpy_bytes += size;
  
  for (i = 1; i < k; i++) {
    galois_region_xor(data_ptrs[i], parity_ptr, size);
    jerasure_total_xor_bytes += size;
  }
}

int jerasure_invert_matrix(int *mat, int *inv, int rows, int w)
{
  int cols, i, j, k, x, rs2;
  int row_start, tmp, inverse;
 
  cols = rows;

  k = 0;
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      inv[k] = (i == j) ? 1 : 0;
      k++;
    }
  }

  /* First -- convert into upper triangular  */
  for (i = 0; i < cols; i++) {
    row_start = cols*i;

    /* Swap rows if we ave a zero i,i element.  If we can't swap, then the 
       matrix was not invertible  */

    if (mat[row_start+i] == 0) { 
      for (j = i+1; j < rows && mat[cols*j+i] == 0; j++) ;
      if (j == rows) return -1;
      rs2 = j*cols;
      for (k = 0; k < cols; k++) {
        tmp = mat[row_start+k];
        mat[row_start+k] = mat[rs2+k];
        mat[rs2+k] = tmp;
        tmp = inv[row_start+k];
        inv[row_start+k] = inv[rs2+k];
        inv[rs2+k] = tmp;
      }
    }
 
    /* Multiply the row by 1/element i,i  */
    tmp = mat[row_start+i];
    if (tmp != 1) {
      inverse = galois_single_divide(1, tmp, w);
      for (j = 0; j < cols; j++) { 
        mat[row_start+j] = galois_single_multiply(mat[row_start+j], inverse, w);
        inv[row_start+j] = galois_single_multiply(inv[row_start+j], inverse, w);
      }
    }

    /* Now for each j>i, add A_ji*Ai to Aj  */
    k = row_start+i;
    for (j = i+1; j != cols; j++) {
      k += cols;
      if (mat[k] != 0) {
        if (mat[k] == 1) {
          rs2 = cols*j;
          for (x = 0; x < cols; x++) {
            mat[rs2+x] ^= mat[row_start+x];
            inv[rs2+x] ^= inv[row_start+x];
          }
        } else {
          tmp = mat[k];
          rs2 = cols*j;
          for (x = 0; x < cols; x++) {
            mat[rs2+x] ^= galois_single_multiply(tmp, mat[row_start+x], w);
            inv[rs2+x] ^= galois_single_multiply(tmp, inv[row_start+x], w);
          }
        }
      }
    }
  }

  /* Now the matrix is upper triangular.  Start at the top and multiply down  */

  for (i = rows-1; i >= 0; i--) {
    row_start = i*cols;
    for (j = 0; j < i; j++) {
      rs2 = j*cols;
      if (mat[rs2+i] != 0) {
        tmp = mat[rs2+i];
        mat[rs2+i] = 0; 
        for (k = 0; k < cols; k++) {
          inv[rs2+k] ^= galois_single_multiply(tmp, inv[row_start+k], w);
        }
      }
    }
  }
  return 0;
}

int jerasure_invertible_matrix(int *mat, int rows, int w)
{
  int cols, i, j, k, x, rs2;
  int row_start, tmp, inverse;
 
  cols = rows;

  /* First -- convert into upper triangular  */
  for (i = 0; i < cols; i++) {
    row_start = cols*i;

    /* Swap rows if we ave a zero i,i element.  If we can't swap, then the 
       matrix was not invertible  */

    if (mat[row_start+i] == 0) { 
      for (j = i+1; j < rows && mat[cols*j+i] == 0; j++) ;
      if (j == rows) return 0;
      rs2 = j*cols;
      for (k = 0; k < cols; k++) {
        tmp = mat[row_start+k];
        mat[row_start+k] = mat[rs2+k];
        mat[rs2+k] = tmp;
      }
    }
 
    /* Multiply the row by 1/element i,i  */
    tmp = mat[row_start+i];
    if (tmp != 1) {
      inverse = galois_single_divide(1, tmp, w);
      for (j = 0; j < cols; j++) { 
        mat[row_start+j] = galois_single_multiply(mat[row_start+j], inverse, w);
      }
    }

    /* Now for each j>i, add A_ji*Ai to Aj  */
    k = row_start+i;
    for (j = i+1; j != cols; j++) {
      k += cols;
      if (mat[k] != 0) {
        if (mat[k] == 1) {
          rs2 = cols*j;
          for (x = 0; x < cols; x++) {
            mat[rs2+x] ^= mat[row_start+x];
          }
        } else {
          tmp = mat[k];
          rs2 = cols*j;
          for (x = 0; x < cols; x++) {
            mat[rs2+x] ^= galois_single_multiply(tmp, mat[row_start+x], w);
          }
        }
      }
    }
  }
  return 1;
}

/* Converts a list-style version of the erasures into an array of k+m elements
   where the element = 1 if the index has been erased, and zero otherwise */

int *jerasure_erasures_to_erased(int k, int m, int *erasures)
{
  int td;
  int t_non_erased;
  int *erased;
  int i;

  td = k+m;
  erased = talloc(int, td);
  if (erased == NULL) return NULL;
  t_non_erased = td;

  for (i = 0; i < td; i++) erased[i] = 0;

  for (i = 0; erasures[i] != -1; i++) {
    if (erased[erasures[i]] == 0) {
      erased[erasures[i]] = 1;
      t_non_erased--;
      if (t_non_erased < k) {
        free(erased);
        return NULL;
      }
    }
  }
  return erased;
}
  
void jerasure_free_schedule(int **schedule)
{
  int i;

  for (i = 0; schedule[i][0] >= 0; i++) free(schedule[i]);
  free(schedule[i]);
  free(schedule);
}

void jerasure_free_schedule_cache(int k, int m, int ***cache)
{
  int e1, e2;

  if (m != 2) {
    fprintf(stderr, "jerasure_free_schedule_cache(): m must equal 2\n");
    assert(0);
  }

  for (e1 = 0; e1 < k+m; e1++) {
    for (e2 = 0; e2 < e1; e2++) {
      jerasure_free_schedule(cache[e1*(k+m)+e2]);
    }
    jerasure_free_schedule(cache[e1*(k+m)+e1]);
  }
  free(cache);
}

void jerasure_matrix_dotprod(int k, int w, int *matrix_row,
                          int *src_ids, int dest_id,
                          char **data_ptrs, char **coding_ptrs, int size)
{
  int init;
  char *dptr, *sptr;
  int i;

  if (w != 1 && w != 8 && w != 16 && w != 32) {
    fprintf(stderr, "ERROR: jerasure_matrix_dotprod() called and w is not 1, 8, 16 or 32\n");
    assert(0);
  }

  init = 0;

  dptr = (dest_id < k) ? data_ptrs[dest_id] : coding_ptrs[dest_id-k];

  /* First copy or xor any data that does not need to be multiplied by a factor */

  for (i = 0; i < k; i++) {
    if (matrix_row[i] == 1) {
      if (src_ids == NULL) {
        sptr = data_ptrs[i];
      } else if (src_ids[i] < k) {
        sptr = data_ptrs[src_ids[i]];
      } else {
        sptr = coding_ptrs[src_ids[i]-k];
      }
      if (init == 0) {
        memcpy(dptr, sptr, size);
        jerasure_total_memcpy_bytes += size;
        init = 1;
      } else {
        galois_region_xor(sptr, dptr, size);
        jerasure_total_xor_bytes += size;
      }
    }
  }

  /* Now do the data that needs to be multiplied by a factor */

  for (i = 0; i < k; i++) {
    if (matrix_row[i] != 0 && matrix_row[i] != 1) {
      if (src_ids == NULL) {
        sptr = data_ptrs[i];
      } else if (src_ids[i] < k) {
        sptr = data_ptrs[src_ids[i]];
      } else {
        sptr = coding_ptrs[src_ids[i]-k];
      }
      switch (w) {
        case 8:  galois_w08_region_multiply(sptr, matrix_row[i], size, dptr, init); break;
        case 16: galois_w16_region_multiply(sptr, matrix_row[i], size, dptr, init); break;
        case 32: galois_w32_region_multiply(sptr, matrix_row[i], size, dptr, init); break;
      }
      jerasure_total_gf_bytes += size;
      init = 1;
    }
  }
}


int jerasure_bitmatrix_decode(int k, int m, int w, int *bitmatrix, int row_k_ones, int *erasures,
                            char **data_ptrs, char **coding_ptrs, int size, int packetsize)
{
  int i;
  int *erased;
  int *decoding_matrix;
  int *dm_ids;
  int edd, *tmpids, lastdrive;
  
  erased = jerasure_erasures_to_erased(k, m, erasures);
  if (erased == NULL) return -1;

  /* See jerasure_matrix_decode for the logic of this routine.  This one works just like
     it, but calls the bitmatrix ops instead */

  lastdrive = k;
    
  edd = 0;
  for (i = 0; i < k; i++) {
    if (erased[i]) {
      edd++;
      lastdrive = i;
    } 
  }

  if (row_k_ones != 1 || erased[k]) lastdrive = k;
  
  dm_ids = NULL;
  decoding_matrix = NULL;
  
  if (edd > 1 || (edd > 0 && (row_k_ones != 1 || erased[k]))) {

    dm_ids = talloc(int, k);
    if (dm_ids == NULL) {
      free(erased);
      return -1;
    }
  
    decoding_matrix = talloc(int, k*k*w*w);
    if (decoding_matrix == NULL) {
      free(erased);
      free(dm_ids);
      return -1;
    }
  
    if (jerasure_make_decoding_bitmatrix(k, m, w, bitmatrix, erased, decoding_matrix, dm_ids) < 0) {
      free(erased);
      free(dm_ids);
      free(decoding_matrix);
      return -1;
    }
  }

  for (i = 0; edd > 0 && i < lastdrive; i++) {
    if (erased[i]) {
      jerasure_bitmatrix_dotprod(k, w, decoding_matrix+i*k*w*w, dm_ids, i, data_ptrs, coding_ptrs, size, packetsize);
      edd--;
    }
  }

  if (edd > 0) {
    tmpids = talloc(int, k);
    if (!tmpids) {
      free(erased);
      free(dm_ids);
      free(decoding_matrix);
      return -1;
    }
    for (i = 0; i < k; i++) {
      tmpids[i] = (i < lastdrive) ? i : i+1;
    }
    jerasure_bitmatrix_dotprod(k, w, bitmatrix, tmpids, lastdrive, data_ptrs, coding_ptrs, size, packetsize);
    free(tmpids);
  }

  for (i = 0; i < m; i++) {
    if (erased[k+i]) {
      jerasure_bitmatrix_dotprod(k, w, bitmatrix+i*k*w*w, NULL, k+i, data_ptrs, coding_ptrs, size, packetsize);
    }
  }

  free(erased);
  if (dm_ids != NULL) free(dm_ids);
  if (decoding_matrix != NULL) free(decoding_matrix);

  return 0;
}

static char **set_up_ptrs_for_scheduled_decoding(int k, int m, int *erasures, char **data_ptrs, char **coding_ptrs)
{
  int ddf, cdf;
  int *erased;
  char **ptrs;
  int i, j, x;

  ddf = 0;
  cdf = 0;
  for (i = 0; erasures[i] != -1; i++) {
    if (erasures[i] < k) ddf++; else cdf++;
  }
  
  erased = jerasure_erasures_to_erased(k, m, erasures);
  if (erased == NULL) return NULL;

  /* Set up ptrs.  It will be as follows:

       - If data drive i has not failed, then ptrs[i] = data_ptrs[i].
       - If data drive i has failed, then ptrs[i] = coding_ptrs[j], where j is the 
            lowest unused non-failed coding drive.
       - Elements k to k+ddf-1 are data_ptrs[] of the failed data drives.
       - Elements k+ddf to k+ddf+cdf-1 are coding_ptrs[] of the failed data drives.

       The array row_ids contains the ids of ptrs.
       The array ind_to_row_ids contains the row_id of drive i.
  
       However, we're going to set row_ids and ind_to_row in a different procedure.
   */
         
  ptrs = talloc(char *, k+m);
  if (!ptrs) {
    free(erased);
    return NULL;
  }

  j = k;
  x = k;
  for (i = 0; i < k; i++) {
    if (erased[i] == 0) {
      ptrs[i] = data_ptrs[i];
    } else {
      while (erased[j]) j++;
      ptrs[i] = coding_ptrs[j-k];
      j++;
      ptrs[x] = data_ptrs[i];
      x++;
    }
  }
  for (i = k; i < k+m; i++) {
    if (erased[i]) {
      ptrs[x] = coding_ptrs[i-k];
      x++;
    }
  }
  free(erased);
  return ptrs;
}

static int set_up_ids_for_scheduled_decoding(int k, int m, int *erasures, int *row_ids, int *ind_to_row)
{
  int ddf, cdf;
  int *erased;
  int i, j, x;

  ddf = 0;
  cdf = 0;
  for (i = 0; erasures[i] != -1; i++) {
    if (erasures[i] < k) ddf++; else cdf++;
  }
  
  erased = jerasure_erasures_to_erased(k, m, erasures);
  if (erased == NULL) return -1;

  /* See set_up_ptrs_for_scheduled_decoding for how these are set */

  j = k;
  x = k;
  for (i = 0; i < k; i++) {
    if (erased[i] == 0) {
      row_ids[i] = i;
      ind_to_row[i] = i;
    } else {
      while (erased[j]) j++;
      row_ids[i] = j;
      ind_to_row[j] = i;
      j++;
      row_ids[x] = i;
      ind_to_row[i] = x;
      x++;
    }
  }
  for (i = k; i < k+m; i++) {
    if (erased[i]) {
      row_ids[x] = i;
      ind_to_row[i] = x;
      x++;
    }
  }
  free(erased);
  return 0;
}

static int **jerasure_generate_decoding_schedule(int k, int m, int w, int *bitmatrix, int *erasures, int smart)
{
  int i, j, x, drive, y, index, z;
  int *decoding_matrix, *inverse, *real_decoding_matrix;
  int *ptr;
  int *row_ids;
  int *ind_to_row;
  int ddf, cdf;
  int **schedule;
  int *b1, *b2;
 
 /* First, figure out the number of data drives that have failed, and the
    number of coding drives that have failed: ddf and cdf */

  ddf = 0;
  cdf = 0;
  for (i = 0; erasures[i] != -1; i++) {
    if (erasures[i] < k) ddf++; else cdf++;
  }
  
  row_ids = talloc(int, k+m);
  if (!row_ids) return NULL;
  ind_to_row = talloc(int, k+m);
  if (!ind_to_row) {
    free(row_ids);
    return NULL;
  }

  if (set_up_ids_for_scheduled_decoding(k, m, erasures, row_ids, ind_to_row) < 0) {
    free(row_ids);
    free(ind_to_row);
    return NULL;
  }

  /* Now, we're going to create one decoding matrix which is going to 
     decode everything with one call.  The hope is that the scheduler
     will do a good job.    This matrix has w*e rows, where e is the
     number of erasures (ddf+cdf) */

  real_decoding_matrix = talloc(int, k*w*(cdf+ddf)*w);
  if (!real_decoding_matrix) {
    free(row_ids);
    free(ind_to_row);
    return NULL;
  }

  /* First, if any data drives have failed, then initialize the first
     ddf*w rows of the decoding matrix from the standard decoding
     matrix inversion */

  if (ddf > 0) {
    
    decoding_matrix = talloc(int, k*k*w*w);
    if (!decoding_matrix) {
      free(row_ids);
      free(ind_to_row);
      return NULL;
    }
    ptr = decoding_matrix;
    for (i = 0; i < k; i++) {
      if (row_ids[i] == i) {
        bzero(ptr, k*w*w*sizeof(int));
        for (x = 0; x < w; x++) {
          ptr[x+i*w+x*k*w] = 1;
        } 
      } else {
        memcpy(ptr, bitmatrix+k*w*w*(row_ids[i]-k), k*w*w*sizeof(int));
      }
      ptr += (k*w*w);
    }
    inverse = talloc(int, k*k*w*w);
    if (!inverse) {
      free(row_ids);
      free(ind_to_row);
      free(decoding_matrix);
      return NULL;
    }
    jerasure_invert_bitmatrix(decoding_matrix, inverse, k*w);

/*    printf("\nMatrix to invert\n");
    jerasure_print_bitmatrix(decoding_matrix, k*w, k*w, w);
    printf("\n");
    printf("\nInverse\n");
    jerasure_print_bitmatrix(inverse, k*w, k*w, w);
    printf("\n"); */

    free(decoding_matrix);
    ptr = real_decoding_matrix;
    for (i = 0; i < ddf; i++) {
      memcpy(ptr, inverse+k*w*w*row_ids[k+i], sizeof(int)*k*w*w);
      ptr += (k*w*w);
    }
    free(inverse);
  } 

  /* Next, here comes the hard part.  For each coding node that needs
     to be decoded, you start by putting its rows of the distribution
     matrix into the decoding matrix.  If there were no failed data
     nodes, then you're done.  However, if there have been failed
     data nodes, then you need to modify the columns that correspond
     to the data nodes.  You do that by first zeroing them.  Then
     whereever there is a one in the distribution matrix, you XOR
     in the corresponding row from the failed data node's entry in
     the decoding matrix.  The whole process kind of makes my head
     spin, but it works.
   */

  for (x = 0; x < cdf; x++) {
    drive = row_ids[x+ddf+k]-k;
    ptr = real_decoding_matrix + k*w*w*(ddf+x);
    memcpy(ptr, bitmatrix+drive*k*w*w, sizeof(int)*k*w*w);

    for (i = 0; i < k; i++) {
      if (row_ids[i] != i) {
        for (j = 0; j < w; j++) {
          bzero(ptr+j*k*w+i*w, sizeof(int)*w);
        }
      }  
    }

    /* There's the yucky part */

    index = drive*k*w*w;
    for (i = 0; i < k; i++) {
      if (row_ids[i] != i) {
        b1 = real_decoding_matrix+(ind_to_row[i]-k)*k*w*w;
        for (j = 0; j < w; j++) {
          b2 = ptr + j*k*w;
          for (y = 0; y < w; y++) {
            if (bitmatrix[index+j*k*w+i*w+y]) {
              for (z = 0; z < k*w; z++) {
                b2[z] = b2[z] ^ b1[z+y*k*w];
              }
            }
          }
        }
      }  
    }
  }

/*
  printf("\n\nReal Decoding Matrix\n\n");
  jerasure_print_bitmatrix(real_decoding_matrix, (ddf+cdf)*w, k*w, w);
  printf("\n"); */
  if (smart) {
    schedule = jerasure_smart_bitmatrix_to_schedule(k, ddf+cdf, w, real_decoding_matrix);
  } else {
    schedule = jerasure_dumb_bitmatrix_to_schedule(k, ddf+cdf, w, real_decoding_matrix);
  }
  free(row_ids);
  free(ind_to_row);
  free(real_decoding_matrix);
  return schedule;
}

int jerasure_schedule_decode_lazy(int k, int m, int w, int *bitmatrix, int *erasures,
                            char **data_ptrs, char **coding_ptrs, int size, int packetsize, 
                            int smart)
{
  int i, tdone;
  char **ptrs;
  int **schedule;
 
  ptrs = set_up_ptrs_for_scheduled_decoding(k, m, erasures, data_ptrs, coding_ptrs);
  if (ptrs == NULL) return -1;

  schedule = jerasure_generate_decoding_schedule(k, m, w, bitmatrix, erasures, smart);
  if (schedule == NULL) {
    free(ptrs);
    return -1;
  }

  for (tdone = 0; tdone < size; tdone += packetsize*w) {
  jerasure_do_scheduled_operations(ptrs, schedule, packetsize);
    for (i = 0; i < k+m; i++) ptrs[i] += (packetsize*w);
  }

  jerasure_free_schedule(schedule);
  free(ptrs);

  return 0;
}

int jerasure_schedule_decode_cache(int k, int m, int w, int ***scache, int *erasures,
                            char **data_ptrs, char **coding_ptrs, int size, int packetsize)
{
  int i, tdone;
  char **ptrs;
  int **schedule;
  int index;
 
  if (erasures[1] == -1) {
    index = erasures[0]*(k+m) + erasures[0];
  } else if (erasures[2] == -1) {
    index = erasures[0]*(k+m) + erasures[1];
  } else {
    return -1;
  }

  schedule = scache[index];

  ptrs = set_up_ptrs_for_scheduled_decoding(k, m, erasures, data_ptrs, coding_ptrs);
  if (ptrs == NULL) return -1;


  for (tdone = 0; tdone < size; tdone += packetsize*w) {
  jerasure_do_scheduled_operations(ptrs, schedule, packetsize);
    for (i = 0; i < k+m; i++) ptrs[i] += (packetsize*w);
  }

  free(ptrs);

  return 0;
}

/* This only works when m = 2 */

int ***jerasure_generate_schedule_cache(int k, int m, int w, int *bitmatrix, int smart)
{
  int ***scache;
  int erasures[3];
  int e1, e2;
 
  /* Ok -- this is yucky, but it's how I'm doing it.  You will make an index out
     of erasures, which will be  e1*(k+m)+(e2).  If there is no e2, then e2 = e1.
     Isn't that clever and confusing.  Sorry.

     We're not going to worry about ordering -- in other words, the schedule for
     e1,e2 will be the same as e2,e1.  They will have the same pointer -- the 
     schedule will not be duplicated. */

  if (m != 2) return NULL;

  scache = talloc(int **, (k+m)*(k+m+1));
  if (scache == NULL) return NULL;
  
  for (e1 = 0; e1 < k+m; e1++) {
    erasures[0] = e1;
    for (e2 = 0; e2 < e1; e2++) {
      erasures[1] = e2;
      erasures[2] = -1;
      scache[e1*(k+m)+e2] = jerasure_generate_decoding_schedule(k, m, w, bitmatrix, erasures, smart);
      scache[e2*(k+m)+e1] = scache[e1*(k+m)+e2];
    }
    erasures[1] = -1;
    scache[e1*(k+m)+e1] = jerasure_generate_decoding_schedule(k, m, w, bitmatrix, erasures, smart);
  }
  return scache;

}

int jerasure_invert_bitmatrix(int *mat, int *inv, int rows)
{
  int cols, i, j, k;
  int tmp;
 
  cols = rows;

  k = 0;
  for (i = 0; i < rows; i++) {
    for (j = 0; j < cols; j++) {
      inv[k] = (i == j) ? 1 : 0;
      k++;
    }
  }

  /* First -- convert into upper triangular */

  for (i = 0; i < cols; i++) {

    /* Swap rows if we have a zero i,i element.  If we can't swap, then the 
       matrix was not invertible */

    if ((mat[i*cols+i]) == 0) { 
      for (j = i+1; j < rows && (mat[j*cols+i]) == 0; j++) ;
      if (j == rows) return -1;
      for (k = 0; k < cols; k++) {
        tmp = mat[i*cols+k]; mat[i*cols+k] = mat[j*cols+k]; mat[j*cols+k] = tmp;
        tmp = inv[i*cols+k]; inv[i*cols+k] = inv[j*cols+k]; inv[j*cols+k] = tmp;
      }
    }
 
    /* Now for each j>i, add A_ji*Ai to Aj */
    for (j = i+1; j != rows; j++) {
      if (mat[j*cols+i] != 0) {
        for (k = 0; k < cols; k++) {
          mat[j*cols+k] ^= mat[i*cols+k]; 
          inv[j*cols+k] ^= inv[i*cols+k];
        }
      }
    }
  }

  /* Now the matrix is upper triangular.  Start at the top and multiply down */

  for (i = rows-1; i >= 0; i--) {
    for (j = 0; j < i; j++) {
      if (mat[j*cols+i]) {
        for (k = 0; k < cols; k++) {
          mat[j*cols+k] ^= mat[i*cols+k]; 
          inv[j*cols+k] ^= inv[i*cols+k];
        }
      }
    }
  } 
  return 0;
}

int jerasure_invertible_bitmatrix(int *mat, int rows)
{
  int cols, i, j, k;
  int tmp;
 
  cols = rows;

  /* First -- convert into upper triangular */

  for (i = 0; i < cols; i++) {

    /* Swap rows if we have a zero i,i element.  If we can't swap, then the 
       matrix was not invertible */

    if ((mat[i*cols+i]) == 0) { 
      for (j = i+1; j < rows && (mat[j*cols+i]) == 0; j++) ;
      if (j == rows) return 0;
      for (k = 0; k < cols; k++) {
        tmp = mat[i*cols+k]; mat[i*cols+k] = mat[j*cols+k]; mat[j*cols+k] = tmp;
      }
    }
 
    /* Now for each j>i, add A_ji*Ai to Aj */
    for (j = i+1; j != rows; j++) {
      if (mat[j*cols+i] != 0) {
        for (k = 0; k < cols; k++) {
          mat[j*cols+k] ^= mat[i*cols+k]; 
        }
      }
    }
  }
  return 1;
}

  
int *jerasure_matrix_multiply(int *m1, int *m2, int r1, int c1, int r2, int c2, int w)
{
  int *product, i, j, k;

  product = (int *) malloc(sizeof(int)*r1*c2);
  for (i = 0; i < r1*c2; i++) product[i] = 0;

  for (i = 0; i < r1; i++) {
    for (j = 0; j < c2; j++) {
      for (k = 0; k < r2; k++) {
        product[i*c2+j] ^= galois_single_multiply(m1[i*c1+k], m2[k*c2+j], w);
      }
    }
  }
  return product;
}

void jerasure_get_stats(double *fill_in)
{
  fill_in[0] = jerasure_total_xor_bytes;
  fill_in[1] = jerasure_total_gf_bytes;
  fill_in[2] = jerasure_total_memcpy_bytes;
  jerasure_total_xor_bytes = 0;
  jerasure_total_gf_bytes = 0;
  jerasure_total_memcpy_bytes = 0;
}

void jerasure_do_scheduled_operations(char **ptrs, int **operations, int packetsize)
{
  char *sptr;
  char *dptr;
  int op;

  for (op = 0; operations[op][0] >= 0; op++) {
    sptr = ptrs[operations[op][0]] + operations[op][1]*packetsize;
    dptr = ptrs[operations[op][2]] + operations[op][3]*packetsize;
    if (operations[op][4]) {
/*      printf("%d,%d %d,%d\n", operations[op][0], 
      operations[op][1], 
      operations[op][2], 
      operations[op][3]); 
      printf("xor(0x%x, 0x%x -> 0x%x, %d)\n", sptr, dptr, dptr, packetsize); */
      galois_region_xor(sptr, dptr, packetsize);
      jerasure_total_xor_bytes += packetsize;
    } else {
/*      printf("memcpy(0x%x <- 0x%x)\n", dptr, sptr); */
      memcpy(dptr, sptr, packetsize);
      jerasure_total_memcpy_bytes += packetsize;
    }
  }  
}

void jerasure_schedule_encode(int k, int m, int w, int **schedule,
                                   char **data_ptrs, char **coding_ptrs, int size, int packetsize)
{
  char **ptr_copy;
  int i, tdone;

  ptr_copy = talloc(char *, (k+m));
  for (i = 0; i < k; i++) ptr_copy[i] = data_ptrs[i];
  for (i = 0; i < m; i++) ptr_copy[i+k] = coding_ptrs[i];
  for (tdone = 0; tdone < size; tdone += packetsize*w) {
    jerasure_do_scheduled_operations(ptr_copy, schedule, packetsize);
    for (i = 0; i < k+m; i++) ptr_copy[i] += (packetsize*w);
  }
  free(ptr_copy);
}
    
int **jerasure_dumb_bitmatrix_to_schedule(int k, int m, int w, int *bitmatrix)
{
  int **operations;
  int op;
  int index, optodo, i, j;

  operations = talloc(int *, k*m*w*w+1);
  if (!operations) return NULL;
  op = 0;
  
  index = 0;
  for (i = 0; i < m*w; i++) {
    optodo = 0;
    for (j = 0; j < k*w; j++) {
      if (bitmatrix[index]) {
        operations[op] = talloc(int, 5);
	if (!operations[op]) {
	  // -ENOMEM
          goto error;
        }
        operations[op][4] = optodo;
        operations[op][0] = j/w;
        operations[op][1] = j%w;
        operations[op][2] = k+i/w;
        operations[op][3] = i%w;
        optodo = 1;
        op++;
        
      }
      index++;
    }
  }
  operations[op] = talloc(int, 5);
  if (!operations[op]) {
    // -ENOMEM
    goto error;
  }
  operations[op][0] = -1;
  return operations;

error:
  for (i = 0; i <= op; i++) {
    free(operations[op]);
  }
  free(operations);
  return NULL;
}

int **jerasure_smart_bitmatrix_to_schedule(int k, int m, int w, int *bitmatrix)
{
  int **operations;
  int op;
  int i, j;
  int *diff, *from, *b1, *flink, *blink;
  int *ptr, no, row;
  int optodo;
  int bestrow = 0, bestdiff, top;

/*   printf("Scheduling:\n\n");
  jerasure_print_bitmatrix(bitmatrix, m*w, k*w, w); */

  operations = talloc(int *, k*m*w*w+1);
  if (!operations) return NULL;
  op = 0;
  
  diff = talloc(int, m*w);
  if (!diff) {
    free(operations);
    return NULL;
  }
  from = talloc(int, m*w);
  if (!from) {
    free(operations);
    free(diff);
    return NULL;
  }
  flink = talloc(int, m*w);
  if (!flink) {
    free(operations);
    free(diff);
    free(from);
    return NULL;
  }
  blink = talloc(int, m*w);
  if (!blink) {
    free(operations);
    free(diff);
    free(from);
    free(flink);
    return NULL;
  }

  ptr = bitmatrix;

  bestdiff = k*w+1;
  top = 0;
  for (i = 0; i < m*w; i++) {
    no = 0;
    for (j = 0; j < k*w; j++) {
      no += *ptr;
      ptr++;
    }
    diff[i] = no;
    from[i] = -1;
    flink[i] = i+1;
    blink[i] = i-1;
    if (no < bestdiff) {
      bestdiff = no;
      bestrow = i;
    }
  }

  flink[m*w-1] = -1;
  
  while (top != -1) {
    row = bestrow;
    /* printf("Doing row %d - %d from %d\n", row, diff[row], from[row]);  */

    if (blink[row] == -1) {
      top = flink[row];
      if (top != -1) blink[top] = -1;
    } else {
      flink[blink[row]] = flink[row];
      if (flink[row] != -1) {
        blink[flink[row]] = blink[row];
      }
    }

    ptr = bitmatrix + row*k*w;
    if (from[row] == -1) {
      optodo = 0;
      for (j = 0; j < k*w; j++) {
        if (ptr[j]) {
          operations[op] = talloc(int, 5);
          if (!operations[op]) goto error;
          operations[op][4] = optodo;
          operations[op][0] = j/w;
          operations[op][1] = j%w;
          operations[op][2] = k+row/w;
          operations[op][3] = row%w;
          optodo = 1;
          op++;
        }
      }
    } else {
      operations[op] = talloc(int, 5);
      if (!operations[op]) goto error;
      operations[op][4] = 0;
      operations[op][0] = k+from[row]/w;
      operations[op][1] = from[row]%w;
      operations[op][2] = k+row/w;
      operations[op][3] = row%w;
      op++;
      b1 = bitmatrix + from[row]*k*w;
      for (j = 0; j < k*w; j++) {
        if (ptr[j] ^ b1[j]) {
          operations[op] = talloc(int, 5);
          if (!operations[op]) goto error;
          operations[op][4] = 1;
          operations[op][0] = j/w;
          operations[op][1] = j%w;
          operations[op][2] = k+row/w;
          operations[op][3] = row%w;
          optodo = 1;
          op++;
        }
      }
    }
    bestdiff = k*w+1;
    for (i = top; i != -1; i = flink[i]) {
      no = 1;
      b1 = bitmatrix + i*k*w;
      for (j = 0; j < k*w; j++) no += (ptr[j] ^ b1[j]);
      if (no < diff[i]) {
        from[i] = row;
        diff[i] = no;
      }
      if (diff[i] < bestdiff) {
        bestdiff = diff[i];
        bestrow = i;
      }
    }
  }
  
  operations[op] = talloc(int, 5);
  if (!operations[op]) goto error;
  operations[op][0] = -1;
  free(from);
  free(diff);
  free(blink);
  free(flink);

  return operations;

error:
  for (i = 0; i <= op; i++) {
    free(operations[op]);
  }
  free(operations);
  free(from);
  free(diff);
  free(blink);
  free(flink);
  return NULL;
}

void jerasure_bitmatrix_encode(int k, int m, int w, int *bitmatrix,
                            char **data_ptrs, char **coding_ptrs, int size, int packetsize)
{
  int i;

  if (packetsize%sizeof(long) != 0) {
    fprintf(stderr, "jerasure_bitmatrix_encode - packetsize(%d) %c sizeof(long) != 0\n", packetsize, '%');
    assert(0);
  }
  if (size%(packetsize*w) != 0) {
    fprintf(stderr, "jerasure_bitmatrix_encode - size(%d) %c (packetsize(%d)*w(%d))) != 0\n", 
         size, '%', packetsize, w);
    assert(0);
  }

  for (i = 0; i < m; i++) {
    jerasure_bitmatrix_dotprod(k, w, bitmatrix+i*k*w*w, NULL, k+i, data_ptrs, coding_ptrs, size, packetsize);
  }
}

/*
 * Exported function for use by autoconf to perform quick 
 * spot-check.
 */
int jerasure_autoconf_test()
{
  int x = galois_single_multiply(1, 2, 8);
  if (x != 2) {
    return -1;
  }
  return 0;
}

