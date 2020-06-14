#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <limits.h>
#include <string.h>

#include <name.h>

#if !defined(MPI)
typedef int MPI_Fint;
#else
#define EXA_MPI
#endif

#include <exa.h>

int set_env(char *key_val){
  char key[BUFSIZ],val[BUFSIZ];

  char *token=strtok(key_val,"=");
  if(token) strncpy(key,token,BUFSIZ); else return 1;

  token=strtok(NULL,"=");
  if(token) strncpy(val,token,BUFSIZ); else return 1;

  printf("GPU: key=%s value=%s\n",key,val);
  return setenv(key,val,0);
}

#define read_state FORTRAN_UNPREFIXED(read_state,READ_STATE)
void read_state(MPI_Fint *fcomm){
  exaExternalComm ce;
  exaComm c;

#if defined(MPI)
  ce=MPI_Comm_f2c(*fcomm);
#else
  ce=*fcomm;
#endif
  exaCommCreate(&c,ce);

  exaInt rank=exaCommRank(c),size=exaCommSize(c);

  char exa_dir[BUFSIZ],exa_hmholtz_dir[BUFSIZ];
  if(rank==0){
    char fname[PATH_MAX];
    if(getcwd(fname,sizeof(fname))==NULL)
      if(rank==0){
        fprintf(stderr,"GPU: Error getting current working directory.\n");
        exit(1);
      }
    strcat(fname,"/.state");
    printf("GPU: state file=%s\n",fname);
    FILE *fp=fopen(fname,"r");

    char line[BUFSIZ];
    while(fgets(line,BUFSIZ,fp)!=NULL){
      int len=strlen(line); line[len-1]='\0';
      if(strstr(line,"EXA_DIR"))
        strncpy(exa_dir,line,BUFSIZ);
      else if(strstr(line,"EXA_HMHOLTZ_DIR")!=NULL)
        strncpy(exa_hmholtz_dir,line,BUFSIZ);
    }
  }

  exaCommBcast(c,exa_dir        ,BUFSIZ,exaChar_t,0);
  exaCommBcast(c,exa_hmholtz_dir,BUFSIZ,exaChar_t,0);

  set_env(exa_dir); set_env(exa_hmholtz_dir);
}
