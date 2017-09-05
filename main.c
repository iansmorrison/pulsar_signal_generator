/*
 * main.c
 *
 * Stand-alone test vector generator for PST testing.
 * Takes a pulsar profile from an ASCII file of nbins Stokes parameters.
 *
 * Ian Morrison
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <fcntl.h>
#include <errno.h>
#include <float.h>
#include <complex.h>
#include <math.h>

#include "test_vec_gen.h"



void usage()
{
    fprintf (stderr,
             "test_vec_gen [options] parameter_file output_file\n"
             " -v          be verbose\n\n");
}


int main (int argc, char **argv)
{
    // ASCII signal file generated by signalgen()
    FILE * ascii_out_fid;
    
    // binary output signal file
    FILE * binary_out_fid;
    
    // ASCII input parameter file
    FILE * par_fid;
    
    // ASCII input pulse profile file
    FILE * ascii_pulse_fid;
    
    // binary pulse profile file (temp)
    FILE * binary_pulse_fid;
    
    float stokes[4];
    double double_stokes[4];
    float val[16];
    int i;
    
    // parameters to be fetched from file
    int nbins;      // number of phase bins in pulse profile
    float period;   // pulsar period (s)
    float t0;       // time offset (s) for first sample
    int Nout;       // number of complex samples to be genrated each series
    int nseries;    // number of series (loops) of sample generation
    int format;     // selection of format (0 for 'complextoreal', 1 for 'complextocomplex')
    int shift;      // 1 for an fftshift before the inverse FFT, 0 otherwise
    float f0;       // centre frequency in MHz)
    float f_sample_out; // sampling frequency of output data (MHz)
    float DM;       // dispersion measure to simulate
    int out_type;   // output type signalgen() should generate - set to 0 (ASCII text) for stand-alone code
    
    // flag set in verbose mode
    unsigned int verbose = 0;
    
    int arg = 0;
    
    while ((arg=getopt(argc,argv,"v")) != -1)
    {
        switch (arg)
        {
            case 'v':
                verbose++;
                break;
                
            default:
                usage ();
                return 0;
        } 
    }
    
    if (verbose)
        fprintf(stderr,"\nPST Test Vector Generator V2\n\n");
    
    // check and parse the command line arguments
    if (argc-optind != 2)
    {
        fprintf(stderr, "\nError: parameter and output filenames must be specified on command line\n\n");
        usage();
        exit(EXIT_FAILURE);
    }
    
    char parameter_filename[256];
    strcpy(parameter_filename, argv[optind]);
    
    char output_filename[256];
    strcpy(output_filename, argv[optind+1]);
    
    struct stat buf;
    if (stat (parameter_filename, &buf) < 0)
    {
        fprintf (stderr, "Error: failed to stat parameter file [%s]: %s\n", parameter_filename, strerror(errno));
        exit(EXIT_FAILURE);
    }
    
    size_t filesize = buf.st_size;
    if (verbose)
        fprintf (stderr, "filesize for %s is %zu bytes\n", parameter_filename, filesize);
    
    
    /***************************/
    
    // Retrieve parameters from parameter file
    if (verbose)
        fprintf (stderr, "reading parameters from file...\n");
    
    par_fid = fopen(parameter_filename, "r");
    
    fscanf(par_fid,"nbins=%d period=%f t0=%f Nout=%d nseries=%d format=%d shift=%d f0=%f f_sample_out=%f DM=%f out_type=%d", &nbins, &period, &t0, &Nout, &nseries, &format, &shift, &f0, &f_sample_out, &DM, &out_type);
    
    fclose(par_fid);
    
    // fprintf(stderr,"nbins=%d period=%f t0=%f Nout=%d nseries=%d format=%d shift=%d f0=%f f_sample_out=%f DM=%f out_type=%d\n", nbins, period, t0, Nout, nseries, format, shift, f0, f_sample_out, DM, out_type);
    
    // Check that out_type is 0, which is necessary for stand-alone code generation
    if (out_type != 0)
    {
        fprintf (stderr, "Error: output type from test_vec_gen() must be 0 for ASCII text\n");
        return (0);
    }
    
    // Create binary version of input pulse profile file
    // (because stand-alone Matlab code cannot read ASCII files)
    
    // Open input pulse profile file
    if (verbose)
        fprintf (stderr, "opening input pulse profile file...\n");
    ascii_pulse_fid = fopen("pulse.txt", "r");
    
    // Open output pulse profile file
    if (verbose)
        fprintf (stderr, "opening output pulse profile file...\n");
    binary_pulse_fid = fopen("pulse.bin", "wb");
    
    while(1)
    {
        fscanf(ascii_pulse_fid,"%e%e%e%e", &stokes[0], &stokes[1], &stokes[2], &stokes[3]);
        
        //printf("\n%.15e %.15e %.15e %.15e\n", stokes[0], stokes[1], stokes[2], stokes[3]);
        
        double_stokes[0] = (double)stokes[0];
        double_stokes[1] = (double)stokes[1];
        double_stokes[2] = (double)stokes[2];
        double_stokes[3] = (double)stokes[3];
        
        if(!feof(ascii_pulse_fid))
        {
            fwrite(&double_stokes, sizeof(double), 4, binary_pulse_fid);
        }
        else
        {
            if (verbose)
                fprintf(stderr, "reached end of file...\n");
            break;
        }
    }

    fclose(ascii_pulse_fid);
    fclose(binary_pulse_fid);
    
    // generate test vector file
    
    if (verbose)
        fprintf (stderr, "starting test vector generation...\n");

    test_vec_gen(nbins, (double)period, (double)t0, Nout, nseries, format, shift, (double)f0, (double)f_sample_out, (double)DM, out_type);
    
    if (verbose)
        fprintf (stderr, "test vector generation completed...\n");
    
    // Create binary version of output test vector file
    // (because stand-alone Matlab code cannot write binary files)
    
    // Open ASCII file generated by test_vec_gen()
    if (verbose)
        fprintf (stderr, "opening ASCII test vector file...\n");
    ascii_out_fid = fopen("test_vector.txt", "r");
    
    // Open binary output test vector file
    if (verbose)
        fprintf (stderr, "opening binary output test vector file...\n");
    binary_out_fid = fopen(output_filename, "wb");
    
    while(1)
    {
        fscanf(ascii_out_fid,"%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e", &val[0], &val[1], &val[2], &val[3], &val[4], &val[5], &val[6], &val[7], &val[8], &val[9], &val[10], &val[11], &val[12], &val[13], &val[14], &val[15]);
        
        /* printf("\n%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", val[0], val[1], val[2], val[3], val[4], val[5], val[6], val[7], val[8], val[9], val[10], val[11], val[12], val[13], val[14], val[15]); */
        
        if(!feof(ascii_out_fid))
        {
            fwrite(&val, sizeof(float), 16, binary_out_fid);
        }
        else
        {
            if (verbose)
                fprintf(stderr, "reached end of file...\n");
            break;
        }
    }
    
    fclose(ascii_out_fid);
    fclose(binary_out_fid);
    
    printf("Done!\n");
    
    return(1);

}
