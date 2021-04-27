/*
 *  SCID-TDSE: Simple 1-electron atomic TDSE solver
 *  Copyright (C) 2015-2021 Serguei Patchkovskii, Serguei.Patchkovskii@mbi-berlin.de
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <dx/dx.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <math.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <assert.h>
/*
 * ===========================================================================
 *
 * SCID-TDSE import filter for OpenDX.
 *
 * ===========================================================================
 *
 */
char *scid_dx_id = "$Id: scid-dx.c,v 1.5 2021/04/26 15:43:01 ps Exp $" ;
/*
 * Data and declarations needed by DX 
 */
extern int     DXConnectToServer(char *host, int port);
extern int     DXInputAvailable (int dxfd);
extern Error   DXCallOutboard   (Error (*module)(Object *, Object *), int dxfd);
static int     dx_fd;
/*
 * Mathematical cut-offs
 */
#define SMALL_EXP         (-60.0/* Neglect smaller exp() arguments */)
/*
 * Data structures
 */
static int        verbose            = 1 ;         /* Debugging level                      */
static char       err_buf[1000] ;                  /* Used for constructing error messages */
static char *     scid_name          = NULL ;      /* Name of the external file            */
static time_t     scid_changed       = 0 ;         /* Time of the last change              */
static FILE *     scid_file          = NULL ;      /* Opened file                          */
static int        sd_nradial ;
static int        sd_nspin ;
static int        sd_lmax ;
static int        sd_mmin ;
static int        sd_mmax ;
static int        timestep ;
static double     rm[3][3] ;                       /* Rotation matrix for lab->local      */
static int        wf_size ;                        /* Size of wavefunction fields         */
static double *   r_tab              = NULL ;      /* Radial grid                         */
static double *   l_wfn_re           = NULL ;      /* Left wavefunction and derivatives   */
static double *   l_grad_re          = NULL ;
static double *   l_lap_re           = NULL ;
static double *   l_wfn_im           = NULL ;
static double *   l_grad_im          = NULL ;
static double *   l_lap_im           = NULL ;
static double *   r_wfn_re           = NULL ;      /* Right wavefunction and derivatives  */
static double *   r_grad_re          = NULL ;
static double *   r_lap_re           = NULL ;
static double *   r_wfn_im           = NULL ;
static double *   r_grad_im          = NULL ;
static double *   r_lap_im           = NULL ;
/*
 * Local functions
 */
extern Error   m_SCID(Object *in, Object *out) ;
extern Error   m_FLATTEN(Object *in, Object *out) ;
extern Error   m_NORM1(Object *in, Object *out) ;
extern Error   m_WRAP(Object *in, Object *out) ;
extern Error   m_UNWRAP(Object *in, Object *out) ;
static Error   parse_wavefunction(Object *name) ;
static Error   do_parse_wavefunction(void) ;
static Error   build_wavefunction(Object *mesh, Object *graph) ;
static Error   really_build_wavefunction(int ngrid, const float *grid, float *wfval) ;
static Error   flatten_field(Object *mesh, Object *slab, int axis, int position) ;
static Error   really_flatten_field(const int ndim, const int axis,const int *counts,double weight, float *src_p, float *dst_p) ;
static double  calculate_integration_weight(const int axis, const int ndim, const float *deltas) ;
static double  calculate_norm_weight(const int ndim, const float *deltas) ;
static void    cube_unwrap_ar07(float *p, float *abs, const int count[], const int stride[]) ;
static void    cube_rewrap(float *p, const int count[]) ;
/*
 * Initial entry point - initiate service loop, and wait for instructions
 */
int
main (int argc, char *argv[] )
{
  if (argc<3) {
    fprintf(stderr,"%s: This program is an outboard module for OpenDX.\n", argv[0]) ;
    exit(EXIT_FAILURE) ;
    }
  /*
   * Debug options
   */
  if (verbose>0) {
    freopen("scid-dx.out","w",stdout) ;
    freopen("scid-dx.err","w",stderr) ;
    setlinebuf(stdout) ;
    setlinebuf(stderr) ;
    printf ("Connecting to host %s, socket %d\n", argv[1], atoi(argv[2])) ;
    }
  /*
   * Connect to DX server
   */
  dx_fd  = DXConnectToServer (argv[1], atoi(argv[2]));
  if (dx_fd < 0) {
    fprintf (stderr, "couldn't connect to socket %d on host %s\n", atoi(argv[2]), argv[1]);
    exit (-1);
    }
  /*
   * DX service loop
   */
  for(;;) {
    if (DXInputAvailable (dx_fd) <= 0) break ;
    DXCallOutboard (m_SCID, dx_fd) ;
    }
  /*
   * Server requested termination
   */
  exit(EXIT_SUCCESS) ;
  }

/*
 * This is out module. We take the following inputs and outputs:
 *
 * Inputs:
 *   0 = file name of the SCID-TDSE .data fiel
 *   1 = mesh on which orbitals must be evaluated
 *
 * Outputs:
 *   0 = mesh (same as passed in #3) with left and right wavefunctions attached
 *   1 = comment in the .data file
 */
Error
m_SCID(Object *in, Object *out) 
{
  if (verbose>0) {
    setlinebuf(stdout) ;
    setlinebuf(stderr) ;
    }
  if (verbose>=0) {
    printf("SCID-TDSE import module: %s\n", scid_dx_id ) ;
    fflush(stdout) ;
    }
  /*
   * Input 0 - checkpoint file name. 
   */
  if (verbose>0) printf("Parsing checkpoint file\n") ;
  if (parse_wavefunction(in+0)==ERROR) {
    DXSetError(ERROR_MISSING_DATA,err_buf) ;
    return ERROR ;
    }
  if (verbose>0) printf("Done parsing checkpoint file\n") ;
  /*
   * Input 1 Grid.
   */
  if (in[1]==NULL) {
    DXSetError(ERROR_MISSING_DATA,err_buf) ;
    return ERROR ;
    }
  if (verbose>0) printf("Building orbitals\n") ;
  if (build_wavefunction(in+1,out+0)==ERROR) {
    DXSetError(ERROR_INTERNAL,err_buf) ;
    return ERROR ;
    }
  if (verbose>0) printf("Leaving SCID-TDSE module\n") ;
  return OK ;
  }

static int
wfn_ind(int i, int s, int l, int m)
{
  int ind ;

  if (verbose>=2) {
    if ( (i<=0) || (i>sd_nradial) || (s<=0) || (s>sd_nspin) || (l<0) || (l>sd_lmax) || (m<sd_mmin) || (m>sd_mmax) ) {
      printf("wfn_ind: %d %d %d %d is out of range\n",i,s,l,m) ;
      fflush(stdout) ;
      abort() ;
      }
    }
  ind = (i-1)+((s-1)+(l+(m-sd_mmin)*(sd_lmax+1))*sd_nspin)*sd_nradial ;
  if ( (ind<0) || (ind>wf_size) ) {
    printf("wfn_ind: %d is out of range\n",ind) ;
    fflush(stdout) ;
    abort() ;
    }
  return ind ;
  }

static Error
parse_wavefunction(Object *name)
{
  char        *file ;
  int         reparse = 0 ;
  struct stat stat_buf ;

  if (*name==NULL) {
    strncpy(err_buf,"Missing SCID-TDSE checkpoint filename",sizeof(err_buf)) ;
    return ERROR ;
    }
  if (!DXExtractString(*name, &file)) {
    strncpy(err_buf,"File name is not a string?!",sizeof(err_buf)) ;
    return ERROR ;
    }
  if (scid_name!=NULL) {
    if (strcmp(file,scid_name)!=0) { reparse = 1 ; }
    }
  if (scid_file!=NULL) {
    if (stat(scid_name,&stat_buf)!=0) {
      if (errno==ENOENT) { reparse = 1 ; }
      else {
        snprintf(err_buf,sizeof(err_buf),"fstat: %s: %s",scid_name,strerror(errno)) ;
        return ERROR ;
        }
      }
    if (scid_changed!=stat_buf.st_mtime) { reparse = 1 ; }
    }
  else { reparse = 1 ; }
  if (!reparse) return OK ;
  /*
   *  We only get here if the file must be parsed.
   */
  if (scid_file!=NULL) { fclose(scid_file) ; scid_file = NULL ; }
  if (scid_name!=NULL) { free(scid_name) ; scid_name = NULL ; }
  scid_name = strdup(file) ;
  scid_file = fopen(scid_name,"r") ;
  if (scid_file==NULL) {
    snprintf(err_buf,sizeof(err_buf),"fopen: %s: %s",scid_name,strerror(errno)) ;
    return ERROR ;
    }
  if (fstat(fileno(scid_file),&stat_buf)!=0) {
    snprintf(err_buf,sizeof(err_buf),"fstat: %s: %s",scid_name,strerror(errno)) ;
    return ERROR ;
    }
  scid_changed = stat_buf.st_mtime ;
  /*
   *  Actually parse the file
   */
  return do_parse_wavefunction() ;
  }

static void
free_memory(void)
{
  if(r_tab    !=NULL) { free(r_tab    ) ; r_tab     = NULL ; }
  if(l_wfn_re !=NULL) { free(l_wfn_re ) ; l_wfn_re  = NULL ; }
  if(l_grad_re!=NULL) { free(l_grad_re) ; l_grad_re = NULL ; }
  if(l_lap_re !=NULL) { free(l_lap_re ) ; l_lap_re  = NULL ; }
  if(l_wfn_im !=NULL) { free(l_wfn_im ) ; l_wfn_im  = NULL ; }
  if(l_grad_im!=NULL) { free(l_grad_im) ; l_grad_im = NULL ; }
  if(l_lap_im !=NULL) { free(l_lap_im ) ; l_lap_im  = NULL ; }
  if(r_wfn_re !=NULL) { free(r_wfn_re ) ; r_wfn_re  = NULL ; }
  if(r_grad_re!=NULL) { free(r_grad_re) ; r_grad_re = NULL ; }
  if(r_lap_re !=NULL) { free(r_lap_re ) ; r_lap_re  = NULL ; }
  if(r_wfn_im !=NULL) { free(r_wfn_im ) ; r_wfn_im  = NULL ; }
  if(r_grad_im!=NULL) { free(r_grad_im) ; r_grad_im = NULL ; }
  if(r_lap_im !=NULL) { free(r_lap_im ) ; r_lap_im  = NULL ; }
  }

static void
allocate_memory (void)
{
  wf_size   = sd_nradial*sd_nspin*(1+sd_lmax)*(sd_mmax-sd_mmin+1) ;
  if (verbose>=2) {
    printf("sd_nradial = %d, sd_nspin = %d, sd_lmax = %d, sd_mmin = %d, sd_mmax = %d\n", sd_nradial, sd_nspin, sd_lmax, sd_mmax, sd_mmin) ;
    printf("wf_size    = %d\n",wf_size) ;
    }
  r_tab     = (double *) calloc( sd_nradial+1, sizeof(double) ) ;
  l_wfn_re  = (double *) calloc( wf_size,      sizeof(double) ) ;
  l_grad_re = (double *) calloc( wf_size,      sizeof(double) ) ;
  l_lap_re  = (double *) calloc( wf_size,      sizeof(double) ) ;
  l_wfn_im  = (double *) calloc( wf_size,      sizeof(double) ) ;
  l_grad_im = (double *) calloc( wf_size,      sizeof(double) ) ;
  l_lap_im  = (double *) calloc( wf_size,      sizeof(double) ) ;
  r_wfn_re  = (double *) calloc( wf_size,      sizeof(double) ) ;
  r_grad_re = (double *) calloc( wf_size,      sizeof(double) ) ;
  r_lap_re  = (double *) calloc( wf_size,      sizeof(double) ) ;
  r_wfn_im  = (double *) calloc( wf_size,      sizeof(double) ) ;
  r_grad_im = (double *) calloc( wf_size,      sizeof(double) ) ;
  r_lap_im  = (double *) calloc( wf_size,      sizeof(double) ) ;
  }

static Error
do_parse_wavefunction(void) 
{
  int ipt, lval, mval, sval ;
  int ind ;

  rewind(scid_file) ;
  if (fscanf(scid_file,"%d%d%d%d%d%d",&sd_nradial,&sd_nspin,&sd_lmax,&sd_mmin,&sd_mmax,&timestep)!=6) {
    strncpy(err_buf,"Error reading SCID-TDSE header",sizeof(err_buf)) ;
    return ERROR ;
    }
  if ( (fscanf(scid_file,"%lg%lg%lg",&rm[0][0],&rm[0][1],&rm[0][2])!=3) || 
       (fscanf(scid_file,"%lg%lg%lg",&rm[1][0],&rm[1][1],&rm[1][2])!=3) || 
       (fscanf(scid_file,"%lg%lg%lg",&rm[2][0],&rm[2][1],&rm[2][2])!=3) ) {
    strncpy(err_buf,"Error reading SCID-TDSE rotation matrix",sizeof(err_buf)) ;
    return ERROR ;
    }
  free_memory () ;
  allocate_memory () ;
  if ( (r_tab    == NULL)|| (l_wfn_re == NULL)|| (l_grad_re== NULL)|| (l_lap_re == NULL)||
       (l_wfn_im == NULL)|| (l_grad_im== NULL)|| (l_lap_im == NULL)|| (r_wfn_re == NULL)||
       (r_grad_re== NULL)|| (r_lap_re == NULL)|| (r_wfn_im == NULL)|| (r_grad_im== NULL)||
       (r_lap_im == NULL) ) {
    strncpy(err_buf,"Memory allocation failed",sizeof(err_buf)) ;
    return ERROR ;
    }
  for (ipt=1;ipt<=sd_nradial;ipt++) {
    if (fscanf(scid_file,"%lg",&r_tab[ipt])!=1) {
      snprintf(err_buf,sizeof(err_buf),"Reading radial grid element %d",ipt) ;
      return ERROR ;
      }
    }
  for (mval=sd_mmin;mval<=sd_mmax;mval++) {
    for (lval=0;lval<=sd_lmax;lval++) {
      for (sval=1;sval<=sd_nspin;sval++) {
        /* Left WF */
        for (ipt=1;ipt<=sd_nradial;ipt++) {
          ind = wfn_ind(ipt,sval,lval,mval) ;
          if (fscanf(scid_file,"%lg%lg",l_wfn_re+ind,l_wfn_im+ind)!=2) {
            snprintf(err_buf,sizeof(err_buf),"Reading left wfn index %d (%d,%d,%d,%d)",ind,ipt,sval,lval,mval) ;
            return ERROR ;
            }
          }
        /* Left WF grad */
        for (ipt=1;ipt<=sd_nradial;ipt++) {
          ind = wfn_ind(ipt,sval,lval,mval) ;
          if (fscanf(scid_file,"%lg%lg",l_grad_re+ind,l_grad_im+ind)!=2) {
            snprintf(err_buf,sizeof(err_buf),"Reading left wfn grad index %d (%d,%d,%d,%d)",ind,ipt,sval,lval,mval) ;
            return ERROR ;
            }
          }
        /* Left WF lap */
        for (ipt=1;ipt<=sd_nradial;ipt++) {
          ind = wfn_ind(ipt,sval,lval,mval) ;
          if (fscanf(scid_file,"%lg%lg",l_lap_re+ind,l_lap_im+ind)!=2) {
            snprintf(err_buf,sizeof(err_buf),"Reading left wfn lap index %d (%d,%d,%d,%d)",ind,ipt,sval,lval,mval) ;
            return ERROR ;
            }
          }
        /* Right WF */
        for (ipt=1;ipt<=sd_nradial;ipt++) {
          ind = wfn_ind(ipt,sval,lval,mval) ;
          if (fscanf(scid_file,"%lg%lg",r_wfn_re+ind,r_wfn_im+ind)!=2) {
            snprintf(err_buf,sizeof(err_buf),"Reading right wfn index %d (%d,%d,%d,%d)",ind,ipt,sval,lval,mval) ;
            return ERROR ;
            }
          }
        /* Right WF grad */
        for (ipt=1;ipt<=sd_nradial;ipt++) {
          ind = wfn_ind(ipt,sval,lval,mval) ;
          if (fscanf(scid_file,"%lg%lg",r_grad_re+ind,r_grad_im+ind)!=2) {
            snprintf(err_buf,sizeof(err_buf),"Reading right wfn grad index %d (%d,%d,%d,%d)",ind,ipt,sval,lval,mval) ;
            return ERROR ;
            }
          }
        /* Right WF lap */
        for (ipt=1;ipt<=sd_nradial;ipt++) {
          ind = wfn_ind(ipt,sval,lval,mval) ;
          if (fscanf(scid_file,"%lg%lg",r_lap_re+ind,r_lap_im+ind)!=2) {
            snprintf(err_buf,sizeof(err_buf),"Reading right wfn lap index %d (%d,%d,%d,%d)",ind,ipt,sval,lval,mval) ;
            return ERROR ;
            }
          }
        }
      }
    }
  return OK ;
  }

static Error
build_wavefunction(Object *mesh, Object *graph)
{
  Field      top ;
  Array      gr, data ;
  int        rank, shape ;
  Type       type ;
  Category   cat ;
  float *    grid ;
  float *    wfval ;
  int        ngrid ;
  /*
   *  We'll take the mesh given on input, and insert data fields,
   *  representing molecular orbitals.
   */
  top = (Field)DXCopy(*mesh,COPY_STRUCTURE) ;
  if (top==NULL) {
    snprintf(err_buf,sizeof(err_buf),"build_wavefunction: DXCopy: %s", DXGetErrorMessage()) ;
    return ERROR ;
    }
  /*
   *  Get grid parameters
   */
  gr = (Array)DXGetComponentValue(top,"positions") ;
  if (gr==NULL) {
    snprintf(err_buf,sizeof(err_buf),"build_wavefunction: DXGetComponentValue(positions): %s", DXGetErrorMessage()) ;
    DXDelete((Object)top) ;
    return ERROR ;
    }
  if (DXGetArrayInfo(gr,&ngrid,&type,&cat,&rank,&shape)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"build_wavefunction: DXGetArrayInfo(grid): %s", DXGetErrorMessage()) ;
    DXDelete((Object)top) ;
    return ERROR ;
    }
  if ((type!=TYPE_FLOAT)||(cat!=CATEGORY_REAL)||(rank!=1)||(shape!=3)) {
    snprintf(err_buf,sizeof(err_buf),"build_wavefunction: incompatible positions: type = %d cat = %d rank = %d shape = %d",
            (int)type, (int)cat, (int)rank, (int)shape) ;
    DXDelete((Object)top) ;
    return ERROR ;
    }
  /*
   *  Allocate data field to hold the output, and link it in.
   */
  data = DXNewArray(TYPE_FLOAT,CATEGORY_COMPLEX,1,2) ;
  if (data==NULL) {
    snprintf(err_buf,sizeof(err_buf),"build_wavefunction: DXNewArray: %s", DXGetErrorMessage()) ;
    DXDelete((Object)top) ;
    return ERROR ;
    }
  if (DXSetComponentValue(top,"data",(Object)data)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"build_wavefunction: DXSetComponentValue(data): %s", DXGetErrorMessage()) ;
    DXDelete((Object)top) ; DXDelete((Object)data) ;
    return ERROR ;
    }
  if (DXSetStringAttribute((Object)data,"dep","positions")==NULL) {
    snprintf(err_buf,sizeof(err_buf),"build_wavefunction: DXSetStringAttribute: %s", DXGetErrorMessage()) ;
    DXDelete((Object)top) ;
    return ERROR ;
    }
  if (DXAddArrayData(data,0,ngrid,NULL)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"build_wavefunction: DXAddArrayData: %s", DXGetErrorMessage()) ;
    DXDelete((Object)top) ;
    return ERROR ;
    }
  /*
   *  Get grid positions and wavefunction data
   */
  grid  = DXGetArrayData(gr) ;
  wfval = DXGetArrayData(data) ;
  if ( (grid==NULL) || (wfval==NULL) ) {
    snprintf(err_buf,sizeof(err_buf),"build_wavefunction: DXGetArrayData: %s", DXGetErrorMessage()) ;
    DXDelete((Object)top) ;
    return ERROR ;
    }
  /*
   *  That's it for preparatory work - we can now do the numerical part 
   *  and evaluate the wavefunctions
   */
  if (really_build_wavefunction(ngrid,grid,wfval)==ERROR) {
    DXDelete((Object)top) ;
    return ERROR ;
    }
  /*
   * Finalize filed construction, and bail out.
   */
  if (DXEndField(top)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"build_wavefunction: DXEndField: %s", DXGetErrorMessage()) ;
    DXDelete((Object)top) ;
    return ERROR ;
    }
  *graph = (Object) top ;
  return OK ;
  }

static double *
p_ylm(double *ylm, int l_max, int m_min, int m_max, int m, int l)
{
  int ind ;

  if (verbose>=3) {
    if ( (l<0) || (l>l_max) || (m<m_min) || (m>m_max) ) {
      printf("p_ylm: l_max = %d, m_min = %d, m_max = %d, m = %d, l = %d is not in range\n", l_max, m_min, m_max, m, l) ;
      fflush(stdout) ;
      abort() ;
      }
    }
  ind = 2 * ((m-m_min) + (m_max-m_min+1) * l) ;
  if (verbose>=3) {
    if ( (ind<0) || (ind>=(sd_lmax+1)*(sd_mmax-sd_mmin+1)*2-1) ) {
      printf("p_ylm: ind = %d is not in range\n",ind) ;
      fflush(stdout) ;
      abort() ;
      }
    }
  return ylm + ind ;
  }

static void
complex_mul(const double a[2], const double b[2], double c[2])
{
  c[0] = a[0]*b[0] - a[1]*b[1] ;
  c[1] = a[0]*b[1] + a[1]*b[0] ;
  }

static void
MathAllYLM2(int l_max, int m_min, int m_max, const double dir[3], double *ylm)
{
  const double  y00 = 0.2820947917738781434740397 ;
  int           k, l, m, m_low, m_high, ipow ;
  double *      sqrtn ;
  double *      rsqrn ;
  double *      p ;
  double        sinth, costh ;
  double        r, xymod ;
  double        xpy[2] ;
  double        xmy[2] ;
  double        vlm[2] ;
  double        tmp[2], tmp2[2] ;
  /* */
  if (verbose>=3) {
    printf("Direction = %lf, %lf, %lf. Lmax = %d, Mmin = %d, Mmax = %d\n", dir[0], dir[1], dir[2], l_max, m_min, m_max) ;
    }
  /* */
  sqrtn = (double *) calloc(2*l_max+4,sizeof(double)) ;
  rsqrn = (double *) calloc(2*l_max+4,sizeof(double)) ;
  /* */
  r      = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]) ;
  xpy[0] = dir[0] ; xpy[1] =  dir[1] ;
  xmy[0] = dir[0] ; xmy[1] = -dir[1] ;
  xymod = sqrt(dir[0]*dir[0] + dir[1]*dir[1]) ;
  if (xymod>0) {
    xpy[0] /= xymod ; xpy[1] /= xymod ;
    xmy[0] /= xymod ; xmy[1] /= xymod ;
    }
  else {
    xpy[0] = 1 ; xpy[1] = 0 ;
    xmy[0] = 1 ; xmy[1] = 0 ;
    }
  /* */
  costh = dir[2] / r ;
  sinth = xymod / r ;
  if (costh> 1.) costh =  1. ;
  if (costh<-1.) costh = -1. ;
  if (sinth< 0.) sinth =  0. ;
  if (sinth> 1.) sinth =  1. ;
  /* */
  sqrtn[0] = 0 ;
  for (k=1;k<=2*l_max+1;k++) {
    sqrtn[k] = sqrt(k) ;
    rsqrn[k] = 1. / sqrt(k) ;
    }
  if ( (m_min<=0) && (m_max>=0) ) {
    p = p_ylm(ylm,l_max,m_min,m_max,0,0) ;
    p[0] = y00 ; p[1] = 0. ;
    }
  /* */
  vlm[0] = y00 ; vlm[1] = 0. ;
  for (l=1;l<=l_max;l++) {
    tmp[0] = -sqrtn[2*l+1]*rsqrn[2*l]*xpy[1]*sinth ;
    tmp[1] =  sqrtn[2*l+1]*rsqrn[2*l]*xpy[0]*sinth ;
    complex_mul(vlm,tmp,tmp2) ;
    vlm[0] = tmp2[0] ; vlm[1] = tmp2[1] ;
    if ( (m_min<=l) && (m_max>=l) ){
      p = p_ylm(ylm,l_max,m_min,m_max,l,l) ;
      p[0] = vlm[0] ; p[1] = vlm[1] ;
      }
    }
  /* */
  vlm[0] = y00 ; vlm[1] = 0. ;
  for (l=1;l<=l_max;l++) {
    tmp[0] =  sqrtn[2*l+1]*rsqrn[2*l]*xmy[1]*sinth ;
    tmp[1] = -sqrtn[2*l+1]*rsqrn[2*l]*xmy[0]*sinth ;
    complex_mul(vlm,tmp,tmp2) ;
    vlm[0] = tmp2[0] ; vlm[1] = tmp2[1] ;
    if ( (m_min<=-l) && (m_max>=-l) ){
      p = p_ylm(ylm,l_max,m_min,m_max,-l,l) ;
      p[0] = vlm[0] ; p[1] = vlm[1] ;
      }
    }
  /* */
  for (l=1;l<=l_max;l++) {
    m_low  = -(l-1) ;
    if (m_low<m_min) m_low = m_min ;
    m_high =  (l-1) ;
    if (m_high>m_max) m_high = m_max ;
    /* */
    for (m=m_low;m<=m_high;m++) {
      p      = p_ylm(ylm,l_max,m_min,m_max,m,l-1) ;
      vlm[0] = -p[1] * costh*sqrtn[2*l-1] ;
      vlm[1] =  p[0] * costh*sqrtn[2*l-1] ;
      /* */
      if (l>=abs(m)+2) {
        p       = p_ylm(ylm,l_max,m_min,m_max,m,l-2) ;
        vlm[0] += p[0] * rsqrn[2*l-3]*sqrtn[l-1-m]*sqrtn[l-1+m] ;
        vlm[1] += p[1] * rsqrn[2*l-3]*sqrtn[l-1-m]*sqrtn[l-1+m] ;
        }
      /* */
      p    = p_ylm(ylm,l_max,m_min,m_max,m,l) ;
      p[0] = vlm[0] * sqrtn[2*l+1]*rsqrn[l-m]*rsqrn[l+m] ;
      p[1] = vlm[1] * sqrtn[2*l+1]*rsqrn[l-m]*rsqrn[l+m] ;
      }
    }
  /* */
  for (l=0;l<=l_max;l++) {
    for (m=m_min;m<=m_max;m++) {
      p = p_ylm(ylm,l_max,m_min,m_max,m,l) ;
      tmp[0] = p[0] ; tmp[1] = p[1] ;
      ipow = -2*m-l + 4*l_max ; // Force ipow to be positive
      assert (ipow>=0) ;
      ipow %= 4 ;
           if (ipow==0) { }
      else if (ipow==1) { p[0] = -tmp[1] ; p[1] =  tmp[0] ; }
      else if (ipow==2) { p[0] = -tmp[0] ; p[1] = -tmp[1] ; }
      else if (ipow==3) { p[0] =  tmp[1] ; p[1] = -tmp[0] ; }
      else { 
        printf ("L = %d, M = %d, ipow = %d\n", l, m, ipow) ;
        fflush(stdout) ;
        abort() ; 
        }
      if (verbose>=3) {
        printf("L = %d, M = %d, YLM = %lg, %lg\n", l, m, p[0], p[1]) ;
        }
      }
    }
  /* */
  free (sqrtn) ;
  free (rsqrn) ;
  }

static int
find_nearest_r(double r)
{
  int ilow, ihigh, imid ;
  double d1, d2 ;

  ilow = 1 ; ihigh = sd_nradial ;
  while ( ihigh>ilow+1 ) {
    imid = (ihigh+ilow)/2 ;
    if (r_tab[imid]>=r) { ihigh = imid ; }
    if (r_tab[imid]<=r) { ilow  = imid ; }
    }
  d1 = r - r_tab[ilow] ;
  d2 = r - r_tab[ihigh] ;
  if (d1<0) d1 = -d1 ;
  if (d2<0) d2 = -d2 ;
  if (d1<=d2) return ilow ;
  else        return ihigh ;
  }

static void
build_single_wf(const double lab_xyz[3],double wfn[4])
{
  double r, dr ;
  double *ylm ;
  double *p ;
  int    ipt, l, m, s, ind ;
  double xyz[3] ;
  double wfnr[2], wfnl[2] ;
  double tmp[2], tmp2[2] ;

  /* 
   *  We want wavefunction in the laboratory coordinates, but out expansion is
   *  in local. It's easier and faster to rotate the point where we want it ...
   */
  xyz[0] = rm[0][0]*lab_xyz[0] + rm[0][1]*lab_xyz[1] + rm[0][2]*lab_xyz[2] ;
  xyz[1] = rm[1][0]*lab_xyz[0] + rm[1][1]*lab_xyz[1] + rm[1][2]*lab_xyz[2] ;
  xyz[2] = rm[2][0]*lab_xyz[0] + rm[2][1]*lab_xyz[1] + rm[2][2]*lab_xyz[2] ;
  /* */
  wfn[0] = 0 ; wfn[1] = 0 ; wfn[2] = 0 ; wfn[3] = 0 ;
  r = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]) ;
  if (verbose>=3) {
    printf("Distance to origin = %lg\n",r) ;
    }
  if (r>r_tab[sd_nradial]) {
    if (verbose>=3) {
      printf("R is too far, return zero\n") ;
      }
    return ;
    }
  if (verbose>=3) {
    printf("Allocating %d doubles for harmonics\n",(sd_lmax+1)*(sd_mmax-sd_mmin+1)*2) ;
    }
  ylm = calloc( (sd_lmax+1)*(sd_mmax-sd_mmin+1)*2, sizeof(double) ) ;
  if (verbose>=3) {
    printf("YLM is at %lx\n",(long int)ylm) ;
    }
  MathAllYLM2(sd_lmax,sd_mmin,sd_mmax,xyz,ylm) ; // ???
  ipt = find_nearest_r(r) ;
  dr  = r - r_tab[ipt] ;
  if (verbose>=3) {
    printf("Nearest position = %d at r = %lf, delta = %lf\n", ipt, r_tab[ipt], dr) ;
    }
  for (m=sd_mmin;m<=sd_mmax;m++) {
    for (l=0;l<=sd_lmax;l++) {
      if ( (m<-l) || (m>l) ) continue ;
      p = p_ylm(ylm,sd_lmax,sd_mmin,sd_mmax,m,l) ;
      for (s=1;s<=sd_nspin;s++) {
        ind = wfn_ind(ipt,s,l,m) ;
        wfnl[0] = l_wfn_re[ind] + dr * l_grad_re[ind] + 0.5*dr*dr * l_lap_re[ind] ;
        wfnl[1] = l_wfn_im[ind] + dr * l_grad_im[ind] + 0.5*dr*dr * l_lap_im[ind] ;
        wfnr[0] = r_wfn_re[ind] + dr * r_grad_re[ind] + 0.5*dr*dr * r_lap_re[ind] ;
        wfnr[1] = r_wfn_im[ind] + dr * r_grad_im[ind] + 0.5*dr*dr * r_lap_im[ind] ;
        /* */
        tmp2[0] = p[0] ; tmp2[1] = -p[0] ;
        complex_mul(wfnl,tmp2,tmp) ;
        wfn[0] += tmp[0] ; wfn[1] += tmp[1] ;
        /* */
        complex_mul(wfnr,p,tmp) ;
        wfn[2] += tmp[0] ; wfn[3] += tmp[1] ;
        }
      }
    }
  if (r!=0.) {
    wfn[0] /= r ; wfn[1] /= r ; wfn[2] /= r ; wfn[3] /= r ; 
    }
  free (ylm) ;
  }

static Error   
really_build_wavefunction(int ngrid, const float *grid, float *wfval)
{
  int         ipt ;
  const float *p_g ;
  float       *p_wf ;
  double      loc_xyz[3] ;
  double      loc_wfn[4] ;

  #pragma omp parallel for default(none) shared(ngrid,grid,wfval,verbose) private(ipt,p_g,p_wf,loc_xyz,loc_wfn)
  for (ipt=0;ipt<ngrid;ipt++) {
    p_g  = grid + 3*ipt ;
    p_wf = wfval + 4*ipt ;
    loc_xyz[0] = p_g[0] ; loc_xyz[1] = p_g[1] ; loc_xyz[2] = p_g[2] ; 
    if (verbose>=3) {
      printf("Looking at grid point %lf, %lf, %lf\n", loc_xyz[0], loc_xyz[1], loc_xyz[2]) ;
      }
    build_single_wf(loc_xyz,loc_wfn) ;
    if (verbose>=3) {
      printf(" Left wfn is %lg, %lg\n", loc_wfn[0], loc_wfn[1]) ;
      printf("Right wfn is %lg, %lg\n", loc_wfn[2], loc_wfn[3]) ;
      }
    p_wf[0] = loc_wfn[0] ; p_wf[1] = loc_wfn[1] ;
    p_wf[2] = loc_wfn[2] ; p_wf[3] = loc_wfn[3] ;
    }
  return OK ;
  }
/*
 *  Ad-hoc analysis modules, copied from gamess-dx.c
 */
/*
 * An ad hoc analysis module - flatten a field along specified direction
 *
 * Inputs:
 *   0 = Field to flatten
 *   1 = Axis to flatten along
 *   2 = Position of the residual slab along the flattening axis
 *
 * Outputs:
 *   0 = Flattened field.
 */
Error
m_FLATTEN(Object *in, Object *out) 
{
  int     axis, position ;
  /*
   *
   */
  if (verbose>0) {
    setlinebuf(stdout) ;
    setlinebuf(stderr) ;
    }
  /*
   * Input 1 - Flattening axis
   */
  if (in[1]==NULL) {
    axis = 0 ;
    }
  else {
    if (!DXExtractInteger(in[1],&axis)) {
      snprintf(err_buf,sizeof(err_buf),"FLATTEN: Error decoding axis argument: %s",DXGetErrorMessage()) ;
      DXSetError(ERROR_MISSING_DATA,err_buf) ;
      return ERROR ;
      }
    }
  /*
   * Input 2 - Slab position
   */
  if (in[2]==NULL) {
    position = -1 ;
    }
  else {
    if (!DXExtractInteger(in[2],&position)) {
      snprintf(err_buf,sizeof(err_buf),"FLATTEN: Error decoding slab argument: %s",DXGetErrorMessage()) ;
      DXSetError(ERROR_MISSING_DATA,err_buf) ;
      return ERROR ;
      }
    }
  /*
   *  Input 0 and output 0 - field to be flattened
   */
  if (flatten_field(in+0,out+0,axis,position)==ERROR) {
    DXSetError(ERROR_INTERNAL,err_buf) ;
    return ERROR ;
    }
  return OK ;
  }
/*
 *  Flatten the field by integrating it
 */
#define MAX_DIMENSIONS (3)
static Error   
flatten_field(Object *mesh, Object *slab, int axis, int position)
{
  Field         top ;
  Array         src_gr, dst_gr, dst_mesh, src_data, dst_data ;
  int           ndim, isrc, idst, ic, dst_npts ;
  int           src_counts[MAX_DIMENSIONS] ;
  float         src_origin[MAX_DIMENSIONS] ;
  int           dst_counts[MAX_DIMENSIONS] ;
  float         dst_origin[MAX_DIMENSIONS] ;
  float         deltas    [MAX_DIMENSIONS*MAX_DIMENSIONS] ;
  double        weight ;
  float *       src_p ;
  float *       dst_p ;
  /*
   *  We'll take the mesh given on input, and insert data fields,
   *  representing molecular orbitals.
   */
  top = (Field)DXCopy(*mesh,COPY_STRUCTURE) ;
  if (top==NULL) {
    snprintf(err_buf,sizeof(err_buf),"flatten_field: DXCopy: %s", DXGetErrorMessage()) ;
    return ERROR ;
    }
  *slab = (Object) top ;
  /*
   *  Get source grid parameters
   */
  src_gr   = (Array)DXGetComponentValue(top,"positions") ;
  src_data = (Array)DXGetComponentValue(top,"data") ;
  if ((src_gr==NULL)||(src_data==NULL)) {
    snprintf(err_buf,sizeof(err_buf),"flatten_field: DXGetComponentValue(positions): %s", 
             DXGetErrorMessage()) ;
    return ERROR ;
    }
  /*
   *  Make sure source grid is of the type we can handle.
   */
  if (DXQueryGridPositions(src_gr,&ndim,NULL,NULL,NULL)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"flatten_field: DXQueryGridPositions: %s", DXGetErrorMessage()) ;
    return ERROR ;
    }
  if (verbose>=1) {
    printf("Source grid contains %d dimensions\n",ndim) ;
    }
  if (ndim<=0||ndim>MAX_DIMENSIONS) {
    snprintf(err_buf,sizeof(err_buf),"flatten_field: input field dimensions out of range: %d", ndim) ;
    return ERROR ;
    }
  if (axis<0||axis>=ndim) {
    snprintf(err_buf,sizeof(err_buf),"flatten_field: axis choice is out of range: %d", axis) ;
    return ERROR ;
    }
  if (DXQueryGridPositions(src_gr,&ndim,src_counts,src_origin,deltas)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"flatten_field: DXQueryGridPositions: %s", DXGetErrorMessage()) ;
    return ERROR ;
    }
  /*
   *  Report some information about the source grid, and copy parameters to the
   *  output slab data
   */
  if (verbose>=1) {
    printf("Source grid origin is: %f %f %f\n",src_origin[0],src_origin[1],src_origin[2]) ;
    }
  for (isrc=0;isrc<ndim;isrc++) {
    if (verbose>=1) {
      printf("Dimension %d: %d points delta %f %f %f\n",isrc,src_counts[isrc],
                 deltas[isrc*ndim+0],deltas[isrc*ndim+1],deltas[isrc*ndim+2]) ;
      }
    dst_counts[isrc] = src_counts[isrc] ;
    dst_origin[isrc] = src_origin[isrc] ;
    }
  /*
   *  Prepare grid for the output slab - we'll squash the chosen dimension, and
   *  adjust the origin to position the slab where the user requested.
   */
  dst_counts[axis] = 1 ;
  if ((position<0)||(position>=src_counts[axis])) {
    position = src_counts[axis]/2 ;
    }
  for (ic=0;ic<ndim;ic++) {
    dst_origin[ic] += deltas[ic]*position ;
    }
  /*
   *  Build regular mesh for the output
   */
  dst_gr   = DXMakeGridPositionsV  (ndim,dst_counts,dst_origin,deltas) ;
  dst_mesh = DXMakeGridConnectionsV(ndim,dst_counts) ;
  if ((dst_gr==NULL) || (dst_mesh==NULL)) {
    snprintf(err_buf,sizeof(err_buf),"flatten_field: MakeGridPositionsV/MakeGridConnectionsV: %s", DXGetErrorMessage()) ;
    DXDelete((Object)dst_gr) ; DXDelete((Object)dst_mesh) ;
    return ERROR ;
    }
  if ((DXSetComponentValue(top,"positions",  (Object)dst_gr  )==NULL) ||
      (DXSetComponentValue(top,"connections",(Object)dst_mesh)==NULL) ||
      (DXSetComponentValue(top,"box",        NULL            )==NULL) ||
      (DXSetComponentValue(top,"colors",     NULL            )==NULL)) {
    snprintf(err_buf,sizeof(err_buf),"flatten_field: DXSetComponentValue: %s", DXGetErrorMessage()) ;
    DXDelete((Object)dst_gr) ;
    return ERROR ;
    }
  if (DXChangedComponentValues(top,"data")==NULL) {
    snprintf(err_buf,sizeof(err_buf),"flatten_field: DXChangedComponentValues: %s", DXGetErrorMessage()) ;
    return ERROR ;
    }
  /*
   *  Allocate data field to hold the output, and link it in.
   */
  dst_npts = 1 ;
  for (idst=0;idst<ndim;idst++) {
    dst_npts *= dst_counts[idst] ;
    }
  dst_data = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,0) ;
  if (dst_data==NULL) {
    snprintf(err_buf,sizeof(err_buf),"flatten_field: DXNewArray: %s", DXGetErrorMessage()) ;
    return ERROR ;
    }
  if (DXSetComponentValue(top,"data",(Object)dst_data)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"flatten_field: DXSetComponentValue(data): %s", DXGetErrorMessage()) ;
    DXDelete((Object)dst_data) ;
    return ERROR ;
    }
  if (DXAddArrayData(dst_data,0,dst_npts,NULL)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"flatten_field: DXAddArrayData: %s", DXGetErrorMessage()) ;
    return ERROR ;
    }
  /*
   *  Figure out integration weight
   */
  weight = calculate_integration_weight(axis,ndim,deltas) ;
  if (verbose>=1) {
    printf("Integration weight is %lf\n", weight) ;
    }
  /*
   *  Get data fields, and proceed with the reduction
   */
  src_p = DXGetArrayData(src_data) ;
  dst_p = DXGetArrayData(dst_data) ;
  if ( (src_p==NULL)||(dst_p==NULL) ) {
    snprintf(err_buf,sizeof(err_buf),"flatten_field: DXGetArrayData: %s", DXGetErrorMessage()) ;
    return ERROR ;
    }
  really_flatten_field(ndim,axis,src_counts,weight,src_p,dst_p) ;
  /*
   * Finalize filed construction, and bail out.
   */
  if (DXEndField(top)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"flatten_field: DXEndField: %s", DXGetErrorMessage()) ;
    return ERROR ;
    }
  return OK ;
  }

static double  
calculate_integration_weight(const int axis, const int ndim, const float *deltas)
{
  int     idim, ic ;
  double  ab, bb, wgt ;
  double  pw[ndim] ; /* Integration direction */
  /*
   *  We'll construct the vector along integration direction, which is
   *  orthogonal to all other directions. It's length is the integration
   *  weight.
   */
  for (ic=0;ic<ndim;ic++) {
    pw[ic] = deltas[axis*ndim+ic] ;
    }
  for (idim=0;idim<ndim;idim++) {
    if (idim==axis) continue ;
    /*
     * Take scalar products
     */
    ab = bb = 0 ;
    for (ic=0;ic<ndim;ic++) {
      ab += deltas[idim*ndim+ic]*pw[ic] ;
      bb += deltas[idim*ndim+ic]*deltas[idim*ndim+ic] ;
      }
    if (bb==0) continue ;
    /*
     * Subtract
     */
    ab /= bb ;
    for (ic=0;ic<ndim;ic++) {
      pw[ic] -= ab*deltas[idim*ndim+ic] ;
      }
    }
  /*
   * Calculate norm of the linearly independent integration direction
   */
  wgt = 0 ;
  for (ic=0;ic<ndim;ic++) {
    wgt += pw[ic]*pw[ic] ;
    }
  return sqrt(wgt) ;
  }
/*
 *
 */
static Error
really_flatten_field(const int ndim, const int axis, const int *counts,double weight, float *src_p, float *dst_p)
{
  int     ipt[ndim] ;
  int     src_base, dst_base, stride, idim, ip ;
  float * p ;
  double  accu ;
  /*
   *  Calculate stride of the elements along the reduction direction, and
   *  initialize the compound index at the same time.
   */
  stride = 1 ;
  for (idim=0;idim<ndim;idim++) {
    ipt[idim] = 0 ;
    if (idim>axis) stride *= counts[idim] ;
    }
  if (verbose>0) printf("Reduction stride = %d\n",stride) ;
  /*
   *  Loop over all surviving dimensions
   */
  for (;;) {
    /*
     *  Calculate the source and destination base indices
     */
    src_base = 0 ; dst_base = 0 ;
    for (idim=0;idim<ndim;idim++) {
                      src_base = counts[idim]*src_base + ipt[idim] ;
      if (idim!=axis) dst_base = counts[idim]*dst_base + ipt[idim] ;
      }
    if (verbose>3) {
      printf("Current point is: %d %d %d. src base = %d, dst base = %d\n",
             ipt[0],ipt[1],ipt[2],src_base,dst_base) ;
      }
    /*
     *  Do the reduction
     */
    accu = 0 ; p = src_p + src_base ;
    for (ip=counts[axis];ip>0;ip--) {
      accu += *p ;
      p    += stride ;
      }
    dst_p[dst_base] = weight * accu ;
    /*
     *  Increment the compound index - the last index runs fastest
     */
    for (idim=ndim-1;idim>=0;idim--) {
      if (idim==axis) continue ;
      if (++ipt[idim]<counts[idim]) break ;
      ipt[idim] = 0 ;
      if ( (idim==0) || ((idim==1)&&(axis==0)) ) {
        return OK ;
        }
      }
    }
  }
/*
 * An ad hoc analysis module - calculate 1-norm of a regular field
 *
 * Inputs:
 *   0 = Field to calculate norm for
 * Outputs:
 *   0 = The norm, what else?
 */
Error
m_NORM1(Object *mesh, Object *norm)
{
  /*
   *  This module is a knock-off of FLATTEN, except that it is much simpler
   *  - all we need to produce is one scalar.
   */
  Array         src_gr, src_data ;
  int           ndim, isrc, ic, npts ;
  int           src_counts[MAX_DIMENSIONS] ;
  float         src_origin[MAX_DIMENSIONS] ;
  float         deltas    [MAX_DIMENSIONS*MAX_DIMENSIONS] ;
  double        weight, sum ;
  float *       src_p ;
  /*
   *  Get source grid parameters
   */
  src_gr   = (Array)DXGetComponentValue((Field)*mesh,"positions") ;
  src_data = (Array)DXGetComponentValue((Field)*mesh,"data") ;
  if ((src_gr==NULL)||(src_data==NULL)) {
    snprintf(err_buf,sizeof(err_buf),"m_NORM1: DXGetComponentValue(positions): %s", 
             DXGetErrorMessage()) ;
    return ERROR ;
    }
  /*
   *  Make sure source grid is of the type we can handle.
   */
  if (DXQueryGridPositions(src_gr,&ndim,NULL,NULL,NULL)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_NORM1: DXQueryGridPositions: %s", DXGetErrorMessage()) ;
    return ERROR ;
    }
  if (verbose>=1) {
    printf("Source grid contains %d dimensions\n",ndim) ;
    }
  if (ndim<=0||ndim>MAX_DIMENSIONS) {
    snprintf(err_buf,sizeof(err_buf),"m_NORM1: input field dimensions out of range: %d", ndim) ;
    return ERROR ;
    }
  if (DXQueryGridPositions(src_gr,&ndim,src_counts,src_origin,deltas)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_NORM1: DXQueryGridPositions: %s", DXGetErrorMessage()) ;
    return ERROR ;
    }
  /*
   *  Count total number of points
   */
  for (isrc=0,npts=1;isrc<ndim;isrc++) npts *= src_counts[isrc] ;
  /*
   *  Report some information about the source grid.
   */
  if (verbose>=1) {
    printf("Source grid origin is: %f %f %f\n",src_origin[0],src_origin[1],src_origin[2]) ;
    for (isrc=0;isrc<ndim;isrc++) {
      printf("Dimension %d: %d points delta %f %f %f\n",isrc,src_counts[isrc],
                 deltas[isrc*ndim+0],deltas[isrc*ndim+1],deltas[isrc*ndim+2]) ;
      }
    }
  /*
   *  Figure out integration weight
   */
  weight = calculate_norm_weight(ndim,deltas) ;
  if (verbose>=1) {
    printf("Integration weight is %lf\n", weight) ;
    }
  if (weight<=0) { return ERROR ; }
  /*
   *  Get data fields, and proceed with the reduction
   */
  src_p = DXGetArrayData(src_data) ;
  if ( src_p==NULL ) {
    snprintf(err_buf,sizeof(err_buf),"m_NORM1: DXGetArrayData: %s", DXGetErrorMessage()) ;
    return ERROR ;
    }
  for (ic=0,sum=0;ic<npts;ic++,src_p++) {
    sum += *src_p ;
    }
  sum *= weight ;
  /*
   * Return the result, and we are done.
   */
  *norm = (Object)DXMakeFloat(sum) ;
  return OK ;
  }
/*
 *  Return volume of the grid element
 */
static double  
calculate_norm_weight(const int ndim, const float *d)
{
  switch (ndim) {
    case 1:
      return fabs(d[0]) ;
    case 2:
      return fabs(d[0]*d[3] - d[1]*d[2]) ;
    case 3:
      return fabs( d[0]*d[4]*d[8] + d[1]*d[5]*d[6] + d[2]*d[3]*d[7]
                 - d[2]*d[4]*d[6] - d[1]*d[3]*d[8] - d[0]*d[5]*d[7] ) ;
    default:
      snprintf(err_buf,sizeof(err_buf),"calculate_norm_weight: dimensionality %d not supported", ndim) ;
      return -1 ;
    }
  }
/*
 * An ad hoc analysis module - perform unwrap() operation assuming input data is smoothly-varying phase
 *
 * Inputs:
 *   0 = Field to unwrap
 *   1 = Type of unwrap; currently ignored
 * Outputs:
 *   0 = Unwrapped field
 */
Error
m_UNWRAP(Object *src, Object *dst)
{
  /*
   *  This module is a knock-off of FLATTEN, except that it is much simpler
   *  - all we need to produce is one scalar.
   */
  Field         top ;
  Array         src_gr, src_data, abs_data, dst_data ;
  int           ndim, ic, npts ;
  int		stride[MAX_DIMENSIONS] ;
  int           counts[MAX_DIMENSIONS] ;
  int           fancy, pass ;
  float *       src_p ;
  float *       abs_p ;
  float *       dst_p ;
  /*
   *  Fancyness level
   */
  fancy = 0 ;
  if (src[2]!=NULL) {
    if (!DXExtractInteger(src[2],&fancy)) {
      snprintf(err_buf,sizeof(err_buf),"m_UNWRAP: DXExtractInteger for the third argument: %s", DXGetErrorMessage()) ;
      DXSetError(ERROR_MISSING_DATA,err_buf) ;
      return ERROR ;
      }
    }
  /*
   *  Start by making a copy of the structure
   */
  top = (Field)DXCopy(*src,COPY_STRUCTURE) ;
  if (top==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_UNWRAP: DXCopy: %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_MISSING_DATA,err_buf) ;
    return ERROR ;
    }
  *dst = (Object) top ;
  /*
   *  From now on, we'll be operating on the copy of the input structure
   */
  /*
   *  Get source grid parameters
   */
  src_gr   = (Array)DXGetComponentValue(top,"positions") ;
  src_data = (Array)DXGetComponentValue(top,"data") ;
  if ((src_gr==NULL)||(src_data==NULL)) {
    snprintf(err_buf,sizeof(err_buf),"m_UNWRAP: DXGetComponentValue(positions): %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_MISSING_DATA,err_buf) ;
    return ERROR ;
    }
  /*
   *  Make sure source grid is of the type we can handle.
   */
  if (DXQueryGridPositions(src_gr,&ndim,NULL,NULL,NULL)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_UNWRAP: DXQueryGridPositions: %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_MISSING_DATA,err_buf) ;
    return ERROR ;
    }
  if (verbose>=1) {
    printf("Source grid contains %d dimensions\n",ndim) ;
    }
  if (ndim<=0||ndim>MAX_DIMENSIONS) {
    snprintf(err_buf,sizeof(err_buf),"m_UNWRAP: input field dimensions out of range: %d", ndim) ;
    DXSetError(ERROR_INTERNAL,err_buf) ;
    return ERROR ;
    }
  if (DXQueryGridPositions(src_gr,&ndim,counts,NULL,NULL)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_UNWRAP: DXQueryGridPositions: %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_INTERNAL,err_buf) ;
    return ERROR ;
    }
  stride[ndim-1] = 1 ;
  for (ic=ndim-1;ic>0;ic--) { stride[ic-1] = stride[ic]*counts[ic] ; }
  npts = stride[0]*counts[0] ;
  /*
   *  Create a copy of the data field: it will be modified!
   */
  if (DXChangedComponentValues(top,"data")==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_UNWRAP: DXChangedComponentValues: %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_INTERNAL,err_buf) ;
    return ERROR ;
    }
  dst_data = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,0) ;
  if (dst_data==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_UNWRAP: DXNewArray: %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_INTERNAL,err_buf) ;
    return ERROR ;
    }
  if (DXSetComponentValue(top,"data",(Object)dst_data)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_UNWRAP: DXSetComponentValue(data): %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_INTERNAL,err_buf) ;
    DXDelete((Object)dst_data) ;
    return ERROR ;
    }
  if (DXAddArrayData(dst_data,0,npts,NULL)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"flatten_field: DXAddArrayData: %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_INTERNAL,err_buf) ;
    return ERROR ;
    }
  /*
   *  Get data fields, and proceed with the unwrap
   */
  src_p = DXGetArrayData(src_data) ;
  dst_p = DXGetArrayData(dst_data) ;
  if (src_p==NULL|| dst_p==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_UNWRAP: DXGetArrayData: %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_INTERNAL,err_buf) ;
    return ERROR ;
    }
  memcpy(dst_p,src_p,sizeof(float)*npts) ;
  /*
   *  Get a pointer to the absolute values of the data; we need this to tweak the scores
   *  WARNING: We do no checking on the data types and sizes!!!!!!!!
   */
  abs_data = (Array)DXGetComponentValue((Field)src[1],"data") ;
  if (abs_data==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_UNWRAP: DXGetComponentValue(data): %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_MISSING_DATA,err_buf) ;
    return ERROR ;
    }
  abs_p = DXGetArrayData(abs_data) ;
  if (abs_p==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_UNWRAP: DXGetArrayData: %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_INTERNAL,err_buf) ;
    return ERROR ;
    }
  /*
   *  Unwrap implementation is specific to the dimensionality of the problem:
   */
  switch (ndim) {
    default:
      snprintf(err_buf,sizeof(err_buf),"m_UNWRAP: unwrap is not implemented for ndim=%d",ndim) ;
      DXSetError(ERROR_INTERNAL,err_buf) ;
      return ERROR ;
    case 3: 
      cube_unwrap_ar07(dst_p,abs_p,counts,stride) ; 
      for (pass=0;pass<fancy;pass++) {
        cube_rewrap(dst_p,counts) ; 
        cube_unwrap_ar07(dst_p,abs_p,counts,stride) ; 
        }
      break ;
    }
  /*
   * Finalize filed construction, and bail out.
   */
  if (DXEndField(top)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_UNWRAP: DXEndField: %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_INTERNAL,err_buf) ;
    return ERROR ;
    }
  /*
   * We are done.
   */
  return OK ;
  }
/* */
static double
step_phase(double ref, double val)
{
  double	step ;

  step = remainder(val-ref,2*M_PI) ;
  /*
   *  If we are trying to change phase by -Pi, force it to +Pi
   */
  if (abs(step+M_PI)<=1e-4) { step += 2*M_PI ; }
  return ref + step ;
  }
/*
 *  Unwrapping routine following Abdul-Rahman et al, Appl. Opt. 46, 6623 (2007).
 */
typedef struct {
  int	voxel[2] ; /* Combined indices of voxels connected by this edge */
  float	score ;    /* Score of this edge				*/
  } edge_score ;
/**/
static int
cuar07_edge_compare(const void *p1, const void *p2)
{
  float	s1 = ((edge_score*)p1)->score ;
  float	s2 = ((edge_score*)p2)->score ;

       if ( s1 > s2 ) return -1 ;
  else if ( s1 < s2 ) return +1 ;
  else                return  0 ;
  }
/**/
static void
cuar07_set_minus1(int n, int *p) 
{
  while (n-->0) { *p++ = -1 ; }
  }
/**/
static float
cuar07_fold(float x)
{
  x = fmodf(x,2*M_PI) ;
  while (x>  M_PI) { x -= 2*M_PI ; }
  while (x<=-M_PI) { x += 2*M_PI ; }
  return x ;
  }
/**/
static float
cuar07_sd(float *p, int d)
{
  return cuar07_fold(p[-d]-p[0]) - cuar07_fold(p[0]-p[d]) ;
  }
/**/
static void
cuar07_voxel_scores(float *dst_p, float *abs_p, const int nv[], const int dl[], float *scr_p)
{
  int		i0, i1, i2 ;	/* Indices along each coordinate */
  int		o0, o1, o2 ;	/* Current pixel offset          */
  float         grad, g2 ;

  for (i0=0,o0=0 ;i0<nv[0];i0++,o0+=dl[0]) {
  for (i1=0,o1=o0;i1<nv[1];i1++,o1+=dl[1]) {
  for (i2=0,o2=o1;i2<nv[2];i2++,o2+=dl[2]) {
    if (i0==0 || i1==0 || i2==0 || i0==nv[0]-1 || i1==nv[1] || i2==nv[2]) {
      /* 
       * This is an edge point; we'll simply give these points a zero score.
       * As the result, they'll be the last points we deal with.
       */
      scr_p[o2] = 0 ;
      }
    else {
      /*
       * General interior point, calculate second derivative
       */
      g2   = 0 ;
      grad = cuar07_sd(dst_p+o2,dl[0]) ; g2 += grad*grad ;
      grad = cuar07_sd(dst_p+o2,dl[1]) ; g2 += grad*grad ;
      grad = cuar07_sd(dst_p+o2,dl[2]) ; g2 += grad*grad ;
      scr_p[o2] = 1.0/(sqrtf(g2) + 1e-3) ; /* Make sure our scores remain finite    */
      /*
       *  Ar07 stops here; This works well, except at the sharp edges.
       *  I find that reducing the score of near-zero voxels helps with
       *  the edge detection.
       */
      scr_p[o2] *= fabsf(abs_p[o2]) ;      /* Reduce the score for near-zero voxels */
      }
    } } }
  }
/**/
static void
cuar07_add_edge(const float dst_p[], const float scr_p[], const int v0, const int dl, const int ind, const int mx, edge_score *edge)
{
  float g1, g2, g3, gx ;
  int   vb ;
  /*
   *  Ar07 uses simple average of the voxel scores. This works very well if there
   *  are no sharp edges; for the edges, this choice tends to yield a very jagged
   *  boundary. Let's try to reduce the score of connections where a large change
   *  happens.
   */
  edge->voxel[0] = v0 ;
  edge->voxel[1] = v0+dl ;
/*edge->score    = scr_p[v0] + scr_p[v0+dl] ;*/
  edge->score    = scr_p[v0] * scr_p[v0+dl] ; /* Try to bias things towards similar scores */
  /*
   *  We can't do a directional correction unless we have at least 4 points along this dimension
   */
  if ((ind>=1) && (ind<mx)) {
    vb = v0 - dl ;
    g1 = cuar07_fold(dst_p[vb+0*dl]-dst_p[vb+1*dl]) ;
    g2 = cuar07_fold(dst_p[vb+1*dl]-dst_p[vb+2*dl]) ;
    g3 = cuar07_fold(dst_p[vb+2*dl]-dst_p[vb+3*dl]) ;
    gx = (g1-g2)*(g1-g2) + (g3-g2)*(g3-g2) ;
    edge->score *= 1.0/(sqrtf(gx) + 1.0) ;
    }
  }
/**/
static void
cuar07_edge_scores(const float dst_p[], const int nv[], const int dl[], const float scr_p[], edge_score *edge_p)
{
  int		i0, i1, i2 ;	/* Indices along each coordinate */
  int		o0, o1, o2 ;	/* Current pixel offset          */

  for (i0=0,o0=0 ;i0<nv[0];i0++,o0+=dl[0]) {
  for (i1=0,o1=o0;i1<nv[1];i1++,o1+=dl[1]) {
  for (i2=0,o2=o1;i2<nv[2];i2++,o2+=dl[2]) {
    /*
     *  Pick out "upward"-facing edges
     */
    if (i0<nv[0]-1) cuar07_add_edge(dst_p,scr_p,o2,dl[0],i0,nv[0]-1,edge_p++) ;
    if (i1<nv[1]-1) cuar07_add_edge(dst_p,scr_p,o2,dl[1],i1,nv[1]-1,edge_p++) ;
    if (i2<nv[2]-1) cuar07_add_edge(dst_p,scr_p,o2,dl[2],i2,nv[2]-1,edge_p++) ;
    } } }
  }
/**/
static void
cuar07_create_region(const int v0, const int v1, float dst_p[], int boss_p[], int link_p[], int tail_p[], int size_p[])
{
  assert(boss_p[v0]==-1 && boss_p[v1]==-1) ;

  boss_p[v0] = v0 ;  
  boss_p[v1] = v0 ;  
  link_p[v0] = v1 ;
  link_p[v1] = -1 ;
  tail_p[v0] = v1 ; /* Tail pointer is only maintained for the boss voxel */
  size_p[v0] = 2 ;  /* Initial size of a group is always 2                */
  dst_p[v1] = step_phase(dst_p[v0],dst_p[v1]) ;
  }
/**/
static void
cuar07_add_voxel(const int v0, const int v1, float dst_p[], int boss_p[], int link_p[], int tail_p[], int size_p[])
{
  int		boss ;
  int		vx ;

  assert(boss_p[v0]!=-1 && boss_p[v1]==-1) ;
  boss       = boss_p[v0] ;
  boss_p[v1] = boss_p[v0] ;
  /*
   *  Find the end of the v0 chain, and attach v1 there
   */
  vx            = tail_p[boss] ;
  link_p[vx]    = v1 ;
  link_p[v1]    = -1 ;
  tail_p[boss]  = v1 ;
  size_p[boss] += 1 ;
  /*
   *  Adjust phase
   */
  dst_p[v1] = step_phase(dst_p[v0],dst_p[v1]) ;
  }
/**/
static void
cuar07_merge_groups(const int v0, const int v1, float dst_p[], int boss_p[], int link_p[], int tail_p[], int size_p[])
{
  double	delta ; /* Phase shift between two groups */
  int		boss0 ; /* Group ID 0                     */
  int		boss1 ; /* Group ID 1                     */
  int           vx ;    /* Voxel iterator                 */
 
  assert(boss_p[v0]!=-1 && boss_p[v1]!=-1 && boss_p[v0]!=boss_p[v1]) ;
  delta = step_phase(dst_p[v0],dst_p[v1]) - dst_p[v1];
  boss0 = boss_p[v0] ;
  boss1 = boss_p[v1] ;
  /*
   *  Switch the v1 group to use v0 group ID, and adjust the phase
   */
  for (vx=boss1;vx!=-1;vx=link_p[vx]) {
    boss_p[vx]  = boss0 ;
    dst_p [vx] += delta ;
    }
  /*
   *  Find the end of the v0 chain, and attach v1 group there
   */
  vx             = tail_p[boss0] ;
  link_p[vx]     = boss1 ;
  tail_p[boss0]  = tail_p[boss1] ;
  size_p[boss0] += size_p[boss1] ;
  }
/**/
static void
cuar07_merge_regions(const int v[], float dst_p[], int boss_p[], int link_p[], int tail_p[], int size_p[])
{
  int		v0, v1 ;	/* Voxel indices	*/
  int		b0, b1 ; 	/* Boss indices		*/
  /*
   * We need to do slightly different things depending of whether
   * the two voxels already belong to a region or not; choose the
   * right course of action here.
   */
  v0 = v[0] ; b0 = boss_p[v0] ;
  v1 = v[1] ; b1 = boss_p[v1] ;
  if (b0==-1)
    if (b1==-1)     /* Both voxels are not part of a group */
      cuar07_create_region(v0,v1,dst_p,boss_p,link_p,tail_p,size_p) ;
    else            /* v1 is part of a group; v0 is not    */
      cuar07_add_voxel(v1,v0,dst_p,boss_p,link_p,tail_p,size_p) ;
  else  
    if (b1==-1)     /* v0 is part of a group; v1 is not    */
      cuar07_add_voxel(v0,v1,dst_p,boss_p,link_p,tail_p,size_p) ;
    else  
      if (b0!=b1)   /* Voxels belong to different groups   */
        if (size_p[b0]>=size_p[b1])
             cuar07_merge_groups(v0,v1,dst_p,boss_p,link_p,tail_p,size_p) ;
        else cuar07_merge_groups(v1,v0,dst_p,boss_p,link_p,tail_p,size_p) ;
      else ; /* Both voxels are part of the same group; nothing to be done  */
  }
/**/
static void
cube_unwrap_ar07(float dst_p[], float abs_p[], const int nv[], const int dl[])
{
  float	*	scr_p ; 	/* Array for the point quality scores	*/
  int *		boss_p ;	/* Array for the first-member indices	*/
  int *		link_p ;	/* Array for the next-member indices 	*/
  int *		tail_p ;	/* Array for the last-member indices 	*/
  int *		size_p ;	/* Array for the group sizes         	*/
  edge_score *	edge_p ;	/* List of edges			*/
  int		n_voxel ;	/* Total number of voxels to unwrap	*/
  int		n_edge ;	/* Total number of edges 		*/
  int		i_edge ; 	/* Edge being processed now		*/

  n_voxel = nv[0]*nv[1]*nv[2] ;
  n_edge  = 3*n_voxel - nv[0]*nv[1] - nv[0]*nv[2] - nv[1]*nv[2] ;
  scr_p   = (float*)      calloc(n_voxel,sizeof(float)) ;
  boss_p  = (int*)        calloc(n_voxel,sizeof(int)) ;
  link_p  = (int*)        calloc(n_voxel,sizeof(int)) ;
  tail_p  = (int*)        calloc(n_voxel,sizeof(int)) ;
  size_p  = (int*)        calloc(n_voxel,sizeof(int)) ;
  edge_p  = (edge_score*) calloc(n_edge, sizeof(edge_score)) ;
  if (scr_p==NULL || boss_p==NULL || link_p==NULL || tail_p==NULL || size_p==NULL || edge_p==NULL) {
    fprintf(stderr,"cube_unwrap_ar07: Allocation failed, aborting unwrap\n") ;
    if (scr_p !=NULL) free(scr_p) ;
    if (boss_p!=NULL) free(boss_p) ;
    if (link_p!=NULL) free(link_p) ;
    if (tail_p!=NULL) free(tail_p) ;
    if (size_p!=NULL) free(size_p) ;
    if (edge_p!=NULL) free(edge_p) ;
    return ;
    }
  cuar07_voxel_scores(dst_p,abs_p,nv,dl,scr_p) ;  /* Calculates voxel scores           */
  cuar07_edge_scores (dst_p,nv,dl,scr_p,edge_p) ; /* Prepares the list of edge scores  */
  qsort(edge_p,n_edge,sizeof(edge_score),cuar07_edge_compare) ;
  cuar07_set_minus1(n_voxel,boss_p) ;             /* All voxels are initially unlinked */
  cuar07_set_minus1(n_voxel,link_p) ;             /* All voxels are initially unlinked */
  cuar07_set_minus1(n_voxel,tail_p) ;             /* All voxels are initially unlinked */
  cuar07_set_minus1(n_voxel,size_p) ;             /* All voxels are initially unlinked */
  /* Ready to unwrap now */
  for (i_edge=0;i_edge<n_edge;i_edge++) {
    /*
     * If these are distinct regions, merge them
     */
    // printf("Merging %d and %d, score %f\n",edge_p[i_edge].voxel[0],edge_p[i_edge].voxel[1],edge_p[i_edge].score) ;
    cuar07_merge_regions(edge_p[i_edge].voxel,dst_p,boss_p,link_p,tail_p,size_p) ;
    }
  /* Unwrap complete */
  free(scr_p) ; free(boss_p) ; free(link_p) ; free(tail_p) ; free(size_p) ; free(edge_p) ;
  }
/**/
static void
cube_rewrap(float *p, const int nv[])
{
  int	n ;
  float val ;

  n = nv[0]*nv[1]*nv[2] ;
  while (n-->0) {
    val = remainder(*p,2*M_PI) ;
    if (val<=-M_PI+1e-5) val = M_PI ;
    if (val>M_PI) val = M_PI ;
    *p++ = val ;
    }
  }
/*
 * An ad hoc analysis module - perform wrap() operation assuming input data is smoothly-varying phase
 *
 * Inputs:
 *   0 = Field to unwrap
 * Outputs:
 *   0 = Unwrapped field
 */
Error
m_WRAP(Object *src, Object *dst)
{
  /*
   *  This module is a knock-off of FLATTEN, except that it is much simpler
   *  - all we need to produce is one scalar.
   */
  Field         top ;
  Array         src_data, dst_data ;
  int           npts ;
  float *       src_p ;
  float *       dst_p ;
  int		i ;
  float         val ;
  /*
   *  Start by making a copy of the structure
   */
  top = (Field)DXCopy(*src,COPY_STRUCTURE) ;
  if (top==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_WRAP: DXCopy: %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_MISSING_DATA,err_buf) ;
    return ERROR ;
    }
  *dst = (Object) top ;
  /*
   *  From now on, we'll be operating on the copy of the input structure
   */
  /*
   *  Get source data parameters
   */
  src_data = (Array)DXGetComponentValue(top,"data") ;
  if (src_data==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_WRAP: DXGetComponentValue(data): %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_MISSING_DATA,err_buf) ;
    return ERROR ;
    }
  if (DXGetArrayInfo(src_data,&npts,NULL,NULL,NULL,NULL)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_WRAP: DXGetArrayInfo(data): %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_MISSING_DATA,err_buf) ;
    return ERROR ;
    }
  /*
   *  Create a copy of the data field: it will be modified!
   */
  if (DXChangedComponentValues(top,"data")==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_WRAP: DXChangedComponentValues: %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_INTERNAL,err_buf) ;
    return ERROR ;
    }
  dst_data = DXNewArray(TYPE_FLOAT,CATEGORY_REAL,0) ;
  if (dst_data==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_WRAP: DXNewArray: %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_INTERNAL,err_buf) ;
    return ERROR ;
    }
  if (DXSetComponentValue(top,"data",(Object)dst_data)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_WRAP: DXSetComponentValue(data): %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_INTERNAL,err_buf) ;
    DXDelete((Object)dst_data) ;
    return ERROR ;
    }
  if (DXAddArrayData(dst_data,0,npts,NULL)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"flatten_field: DXAddArrayData: %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_INTERNAL,err_buf) ;
    return ERROR ;
    }
  /*
   *  Get data fields, and proceed with the unwrap
   */
  src_p = DXGetArrayData(src_data) ;
  dst_p = DXGetArrayData(dst_data) ;
  if (src_p==NULL|| dst_p==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_WRAP: DXGetArrayData: %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_INTERNAL,err_buf) ;
    return ERROR ;
    }
  /*
   *  Unwrap implementation is specific to the dimensionality of the problem:
   */
  for (i=0;i<npts;i++) {
    val = remainder(*src_p++,2*M_PI) ;
    if (val<=-M_PI+1e-5) val = M_PI ;
    if (val>M_PI) val = M_PI ;
    *dst_p++ = val ;
    }
  /*
   * Finalize filed construction, and bail out.
   */
  if (DXEndField(top)==NULL) {
    snprintf(err_buf,sizeof(err_buf),"m_WRAP: DXEndField: %s", DXGetErrorMessage()) ;
    DXSetError(ERROR_INTERNAL,err_buf) ;
    return ERROR ;
    }
  /*
   * We are done.
   */
  return OK ;
  }
