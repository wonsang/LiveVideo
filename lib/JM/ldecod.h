#include "contributors.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/timeb.h>

#if defined WIN32
  #include <io.h>
#else
  #include <unistd.h>
#endif
#include <sys/stat.h>
#include <fcntl.h>

#include <assert.h>

#include "global.h"
#include "rtpp.h"
#include "memalloc.h"
#include "mbuffer.h"
#include "leaky_bucket.h"
#include "fmo.h"
#include "annexb.h"
#include "output.h"
#include "cabac.h"

#include "erc_api.h"

#define JM          "10 (FRExt)"
#define VERSION     "10.2"
#define EXT_VERSION "(FRExt)"

#define LOGFILE     "log.dec"
#define DATADECFILE "dataDec.txt"
#define TRACEFILE   "trace_dec.txt"

extern objectBuffer_t *erc_object_list;
extern ercVariables_t *erc_errorVar;
extern ColocatedParams *Co_located;

// I have started to move the inp and img structures into global variables.
// They are declared in the following lines.  Since inp is defined in conio.h
// and cannot be overridden globally, it is defined here as input
//
// Everywhere, input-> and img-> can now be used either globally or with
// the local override through the formal parameter mechanism

extern FILE* bits;
extern StorablePicture* dec_picture;

struct inp_par    *input;       //!< input parameters from input configuration file
struct snr_par    *snr;         //!< statistics
struct img_par    *img;         //!< image parameters

//Added by Wonsang You at FEB 8, 2007
extern FILE* bitss;
//Added-End

void JMDecHelpExit ();
//void Configure(int ac, char *av[]);
//int main(int argc, char **argv);
void init(struct img_par *img);
void init_frext(struct img_par *img);
void init_conf(struct inp_par *inp, char *config_filename);
void report(struct inp_par *inp, struct img_par *img, struct snr_par *snr);
//DataPartition *AllocPartition(int n);
void FreePartition (DataPartition *dp, int n);
void malloc_slice(struct inp_par *inp, struct img_par *img);
void free_slice(struct inp_par *inp, struct img_par *img);
int init_global_buffers();
void free_global_buffers();
void open_bitstream_file();

//Added by Wonsang You at FEB 2, 2007
int last_object_tracking_main();
int first_object_tracking_main();
//Added-End

