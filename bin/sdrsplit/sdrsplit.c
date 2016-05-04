/************************************************************************/
/*  sdrsplit:								*/
/*	Demultiplex block-multiplexed SEED Data Record (pre-MiniSEED)	*/
/*	and MiniSEED records into separate files per channel.		*/
/*									*/
/*	Douglas Neuhauser						*/
/*	Seismological Laboratory					*/
/*	University of California, Berkeley				*/
/*	doug@seismo.berkeley.edu					*/
/*									*/
/************************************************************************/

/*
 * Copyright (c) 1996-2000 The Regents of the University of California.
 * All Rights Reserved.
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for educational, research and non-profit purposes,
 * without fee, and without a written agreement is hereby granted,
 * provided that the above copyright notice, this paragraph and the
 * following three paragraphs appear in all copies.
 * 
 * Permission to incorporate this software into commercial products may
 * be obtained from the Office of Technology Licensing, 2150 Shattuck
 * Avenue, Suite 510, Berkeley, CA  94704.
 * 
 * IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
 * FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,
 * INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND
 * ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF CALIFORNIA HAS BEEN
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
 * PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
 * CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
 * UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
 */

/*
Modification History
------------------------------------------------------------------------
Date	    Ver	    Who What
------------------------------------------------------------------------
2010/07/16  1.4.11  DSN	Updated to merge contiguous segments of data per 
	channel that are separated in input file by other non-contiguous 
	data of the same channel.
2010/03/16  1.4.9   DSN	Fixed eof detection (again).
2010/03/09  1.4.8   DSN Fixed eof detection.
2010/03/01  1.4.7   DSN Updated to handle multiple blksizes in input file.
2008/02/21  1.4.6   DSN Added -u option to create unique files.
2004/12/11  1.4.5   DSN	Use network in comparison.
1997/09/15  1.4.1   DSN Added -D option to split data based on day.
*/


#ifndef lint
static const char sccsid[] = "$Id: sdrsplit.c,v 1.20 2013/09/17 23:47:02 doug Exp $ ";
#endif

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/time.h>
#include <math.h>
#include <fnmatch.h>
   
#include "version.h"

#ifndef	DEFAULT_NET
#define	DEFAULT_NET	"XX"
#endif

#define	QUOTE(x)	#x
#define	STRING(x)	QUOTE(x)

char *syntax[] = {
"%s version " VERSION,
"%s  [-n] [-c] [-N netid] [-P] [-C] [-D] [-g tol] [-G tol] [-s] [-r]",
"    [-i] [-d n] [-u] [-v] [-h] [input_file]",
"    where:",
"	-n	    No output - do not write demultiplexed output channels.",
"	-c	    For each channel, write index consisting of station ",
"		    channel, and network to stdout.",
"	-N netid    specifies the  network ID to be used if no network ID is",
"		    found in the SEED data record header.",
"		    The default network_id is " DEFAULT_NET ".",
"	-P	    Create filenames using the Passcal naming convention of",
"			yy.dd.hh.mm.ss.net.station.channel",
"	-C	    Split discontinuous data from a channel into multiple",
"		    files of continuous data.  If this option is not used",
"		    all data from a given station-channel will be written to",
"		    a single station-channel file.",
"	-D	    Split data from a channel into multiple day files based",
"		    on the start day of each record.",
"	-g tol	    Consider data for a single station-channel to be ",
"		    non-continuous if the gap between any adjacent blocks for",
"		    a station/channel exceeds tol ticks.  A tick is 1/10 msec.",
"	-G tol	    Consider data for a single station-channel to be ",
"		    non-continuous if the total time slew for ",
"		    a station/channel exceeds tol ticks.  A tick is 1/10 msec.",
"	-s	    Split output based on non-consecutive SEED sequence numbers.",
"	-r	    Split output channels base on sample rate as well as",
"		    station & channel name.  The rate appears in the filename.",
"		    Used in emergencies to separate different sample rate",
"		    channels inadvertently recorded with the same SEED",
"		    channel name.",
"	-i	    Output on stdout an index of stations, channels, and times",
"		    found in the input file.",
"	-x stat_chan_list",
"		    Extract only records for the specified stations & channels.",
"		    List is a comma-delimited list of channel ids",
"		    which may contain wildcards.  Each channel id may be",
"		    station, station.channel, station.network.channel, or",
"		    station.network.channel.location",
"		    eg BKS.BH?,ARC.*,CMB.BK.H*",
"		    If the list contains wildcards, it must be quoted.",
"	-d n	    Debug option.  Values (which can be added together) are",
"			4 = DEBUG_BLOCK",
"			8 = DEBUG_HDRS",
"	-u	    Generate unique files for channels that have duplicate",
"		    data starting at the same time.",
"	-v	    Verbose option.  Provide additional information msgs.",
"	-h	    Help - prints syntax message.",
"	input_file  Input file containing block-multiplexed SEED data records.",
"		    If no input file is specified, data is read from stdin.",
"  output:",
"	Output file for each channel (or contiguous timeperiod), named:",
"	    station.network.channel.D.yyyy.ddd.hhmmss",
"	or  station.network.rate.channel.D.yyyy.ddd.hhmmss",
NULL };

/************************************************************************/

#include "qlib2.h"

#include "params.h"
#include "st_info.h"
#include "procs.h"

typedef struct _chanlist {
    char    station[8];
    char    network[4];
    char    channel[4];
    char    location[4];
} CHANLIST;

/***************************** External variables   *********************/

#define	    MAX_BLOCKS	100
#define	    BUFSIZE	(MAX_BLKSIZE * MAX_BLOCKS)
char	    input[BUFSIZE];/* input data buffer.			*/
int	    herrno;		/* errno from header routines.		*/
int	    blksize = MAX_BLKSIZE;
int	    debug_option = 0;
int	    ignore_seq = 0;
int	    split_on_seq = 0;
int	    split_on_rate = 0;
int	    split_on_gap = 0;
int	    split_on_day = 0;
int	    gap_thresh = -1;
int	    slew_thresh = -1;
int	    gen_index = 0;
int	    unique_files = 0;
int	    verbose = 0;
int	    ignore_unknown_channels = 0;
int	    passcal_filename = 0;
int	    seq_no = 0;
int	    no_output = 0;	/* don't generate any output file.	*/
int	    make_index = 0;
int	    date_fmt = JULIAN_FMT;
ST_INFO	    st_head;		/* head ptr to channel struct linked list*/
INT_TIME    start_time;
INT_TIME    end_time;
int	    start_flag, end_flag;
FILE	    *info;
char	    *default_net = DEFAULT_NET;
char	    *cmdname;
int	    xtractc;
CHANLIST    **xtractv;

int different_day (INT_TIME begtime1, INT_TIME begtime2);
int first_b2000_seq (DATA_HDR *hdr);
int last_b2000_seq (DATA_HDR *hdr);
int xtractable (DATA_HDR *hdr);
int suspend_stream 
   (ST_INFO	*phead);	/* ptr to stream struct linked list.	*/
int reopen_stream
   (ST_INFO	*st_p);		/* ptr to stream struct to reopen.	*/

/************************************************************************/
/*  main:   main program.						*/
/************************************************************************/
int main (int argc, char **argv)
{
    /* Variables needed for getopt. */
    extern char	*optarg;
    extern int	optind, opterr;
    int		c;
    
    FILE	*fpin;		/* input FILE pointer.			*/
    DATA_HDR	*hdr;		/* pointer to data header.		*/
    ST_INFO	*st_p;		/* pointer to channel structure.	*/
    int		nread;		/* # characters read in read loop.	*/
    char	*infile;	/* input file name.			*/
    int		type;		/* type of input file (cont, dialup).	*/
    char	*in;		/* ptr to current input block.		*/
    int		nblocks = 0;	/* # of input blocks processed.		*/
    char	*p;
    int		at_eof = 0;
    int		nleft = 0;

    /*	Initialization.							*/
    cmdname = tail(*argv);
    type = UNK_HDR_TYPE;	/* unknown input header type.		*/
    fpin = stdin;		/* assume input to come from stdin.	*/
    info = stdout;
    get_my_wordorder();

    /*	Parse command line options.					*/
    while ( (c = getopt(argc,argv,"hncrsPCDuvid:g:G:N:x:")) != -1)
	switch (c) {
	case '?':
	case 'h':   print_syntax(cmdname,syntax,info); exit(0); break;
	case 'n':   ++no_output; break;
	case 'c':   ++make_index; break;
    	case 's':   ++split_on_seq; break;
	case 'r':   ++split_on_rate; break;
	case 'u':   ++unique_files; break;
	case 'v':   ++verbose; break;
	case 'd':   debug_option = atoi(optarg); break;
	case 'P':   ++passcal_filename; break;
	case 'C':   ++split_on_gap; break;
	case 'D':   ++split_on_day; break;
	case 'g':   gap_thresh = atoi(optarg); break;
	case 'G':   slew_thresh = atoi(optarg); break;
	case 'N':   default_net = optarg; break;
	case 'i':   ++gen_index; break;
	case 'x':   parse_chanlist(optarg); break;
	}

    /*	Skip over all options and their arguments.			*/
    argv = &(argv[optind]);
    argc -= optind;

    if (strlen(default_net) != 2) {
	p = malloc (3 * sizeof(char));
	sprintf (p, "%-2.2s", default_net);
	default_net = p;
    }

    /*	May have input file, or input may come from stdin.    */
    if (argc--) {
	infile = *argv++;
	if ( (fpin = fopen (infile, "r")) == NULL) {
	    fprintf (info, "unable to open input file %s\n", infile);
	    exit(1);
	}
    }
    else {
	fpin = stdin;
    }

    /*	Loop reading and processing data until EOF or error.		*/
    /*	If we may not know the blocksize of the data, read the maximum	*/
    /*	blocksize and then determine the actual blocksize on the fly.	*/
    /*  NOTE:  Every block could be a different size.			*/

    nblocks = 1;
    nleft = 0;
    at_eof = 0;
    while (1) {
	nread = blkread (fpin, input+nleft, BUFSIZE-nleft);
	if (nread < 0) {
	    fprintf (info, "Error: reading input file\n");
	}
        else if (nread == 0) {
	    ++at_eof;
	}
	else {
	    nleft += nread;
	}
	if (nleft <= 0) break;
	if (nleft < MAX_BLKSIZE) {
	    /* Pad buffer with zeros to 2*MAX_BLKSIZE */
	    blksize = nleft;
	    memset (input+nleft, 0, 2*MAX_BLKSIZE-nleft);
	    ++at_eof;
	}
	in = input;
	for (; nleft >= MAX_BLKSIZE || (nleft > 0 && at_eof); nleft-=blksize, in+=blksize, ++nblocks) {
	    hdr = decode_hdr (in, nleft);
	    if (hdr == NULL) {
		if (herrno == 0) {
		    fprintf (info, "Warning: non-data block %d\n",  nblocks);
		    blksize = 512;
		    continue;
		}
		else {
		    fprintf (info, "Error: Unknown header error at block %d,\n",
			     nblocks);
		    exit(1);
		}
	    }
	    blksize = hdr->blksize;
	    if (! is_data_hdr_ind(hdr->record_type)) {
		fprintf (info, "Warning: non-data block %d\n",  nblocks);
		free_data_hdr (hdr);
		if (blksize <= 0) blksize = 512;
		continue;
	    }
	    if (blksize == 0) {
		fprintf (info, "unable to determine input blksize\n");
		exit(1);
	    }
	    /* Add network code if necessary so we can use it in comparisons. */
	    if (strlen(hdr->network_id)==0 || strcmp(hdr->network_id,"  ")==0) {
		strncpy(hdr->network_id, default_net, 2);
		capnstr(((SDR_HDR *)in)->network_id,default_net,2);
	    }
	    if (debug(DEBUG_BLOCK | DEBUG_HDRS)) {
		fprintf (info, "block %d, seq_no = %d\n", 
			 nblocks, hdr->seq_no); 
	    }
	    st_p = find_st(hdr);
	    if (st_p != NULL) ++(st_p->num_blk);
	    write_channel (st_p, in, blksize);
	}
	/* If we still have any unprocessed data left in the buffer,	*/
	/* move it to the beginning of the buffer.			*/
	if (nleft > 0 && input != in) {
	    memmove (input, in, nleft);
	}
    }

    /*	Update and flush header information.				*/
    for (st_p = st_head.next; st_p != NULL; st_p = st_p->next) {
	flush_channel (st_p, HDR_FLUSH);
    }
    if (make_index && st_head.next != NULL) putchar ('\n');
    return (0);
}

/************************************************************************/
/*  blkread:								*/
/*	Read a chunk of the file.					*/
/*  returns:								*/
/*	>0:		# bytes being returned.				*/
/*	0:		EOF						*/
/*	-1:		Error						*/
/************************************************************************/
int blkread 
   (FILE	*in,		/* input FILE.				*/
    char	*block,		/* sdr record.				*/
    int		blksize)	/* # of bytes to read (0 -> MAX_BLKSIZE)*/
{
    int	    nread;
    int	    need = blksize;
    int	    tread = 0;

    if (need == 0) need = MAX_BLKSIZE;
    while (need > 0 && ((nread = fread (block, 1, need, in)) > 0) ) {
	need -= nread;
	block += nread;
	tread += nread;
    }
    if (tread > 0) return (tread);
    if (feof(in)) return (0);
    return (-1);
}

/************************************************************************/
/*  alloc_channel:							*/
/*	Allocate all memory required for a channel header.		*/
/************************************************************************/
ST_INFO *alloc_channel ()
{
    ST_INFO	*p;

    p = (ST_INFO *)malloc(sizeof(ST_INFO));
    memset ((void *)p, 0, sizeof(ST_INFO));
    p->databuffer = (int *) malloc(MAX_BLKSIZE*sizeof(int));
    memset ((void *)p->databuffer, 0, sizeof(MAX_BLKSIZE*sizeof(int)));
    return (p);
}

/************************************************************************/
/*  free_channel:							*/
/*	Free all memory required for a channel header.			*/
/************************************************************************/
void free_channel
   (ST_INFO	*p)		/* ptr to stream structure.		*/
{
    ST_INFO	*prev;

    if (p != NULL) {
	if (p->init_hdr != NULL) free_data_hdr (p->init_hdr);
	if (p->sum_hdr != NULL) free_data_hdr (p->sum_hdr);
	if (p->prev_hdr != NULL) free_data_hdr (p->prev_hdr);
	if (p->cur_hdr != NULL) free_data_hdr (p->cur_hdr);
	if (p->databuffer != NULL) free (p->databuffer);

	/*  Unchain this ST_INFO from the linked list.			*/
	/*  We ALWAYS have a prev structure, since there is a dummy	*/
	/*  head for the list, but not always a next struct.		*/
	p->prev->next = p->next;    /* unlink fwd pointer.		*/
	if (p->next != NULL) p->next->prev = p->prev;
	prev = p->prev;
	free (p);
    }
}

/************************************************************************/
/*  find_st:								*/
/*	Find the channel info block for this data channel.		*/
/*	Hang the header on the appropriate place for that channel.	*/
/*	Perform all necessary initialization for a new channel.		*/
/************************************************************************/
ST_INFO *find_st
   (DATA_HDR	*hdr)	    /* ptr to DATA_HDR for record.		*/
{
    ST_INFO	*p;
    int		break_channel = 0;
    int		seconds, usecs;
    double	dslew;

    for (p = st_head.next; p != (ST_INFO *)NULL; p = p->next) {
	if (! ( (strcmp(p->station,hdr->station_id)==0) && 
 	    (strcmp(p->network,hdr->network_id)==0) &&
 	    (strcmp(p->channel,hdr->channel_id)==0) &&
	    (strcmp(p->location,hdr->location_id)==0) &&
 	    ((! split_on_rate) || 
	     ((p->sample_rate == hdr->sample_rate) &&
	      (p->sample_rate_mult == hdr->sample_rate_mult)))) ) continue;

	/* Check for any reason that we should break this channel.	*/
	break_channel = 0;
	if ( p != (ST_INFO *)NULL ) {
	    if ( split_on_seq && (p->cur_hdr->seq_no + 1 != hdr->seq_no) ) {
		continue;
	    }
	    if (split_on_day && different_day (p->cur_hdr->begtime, hdr->begtime)) {
		continue;
	    }
	    switch (p->type) {
	      case 'D':
		time_interval2 (p->cur_hdr->num_samples, p->cur_hdr->sample_rate, 
				p->cur_hdr->sample_rate_mult, &seconds, &usecs);
		dslew = tdiff (hdr->begtime, add_time(p->cur_hdr->begtime, seconds, usecs));
		dslew = (dslew / USECS_PER_TICK);
		if (split_on_gap && (fabs(dslew) > p->gap_thresh)) {
		    ++break_channel;
		}
		if (split_on_gap && (fabs(dslew+p->total_slew) > p->slew_thresh)) {
		    ++break_channel;
		}
		break;
	      case 'L':
		break;
	      case 'O':
		if (split_on_gap) {
		    int prev_b_seq, cur_b_seq;
		    prev_b_seq = last_b2000_seq (p->cur_hdr);
		    cur_b_seq = first_b2000_seq (hdr);
		    if (prev_b_seq < 0 || cur_b_seq < 0 ||
			NEXT2000(prev_b_seq) != cur_b_seq) {
			++break_channel;
		    }
		}
	      default:
		break;
	    }
	    if (break_channel) {
		continue;
	    }
	    else {
		break;
	    }
	}
    }

    if ( p == (ST_INFO *)NULL ) {
	/* Create a new header. */
	/* Put it at the head of the list.  */
	p = alloc_channel();
	p->next = st_head.next;
	st_head.next = p;
	p->prev = &st_head;
	if (p->next != NULL) p->next->prev = p;

	strcpy (p->station, hdr->station_id);
	strncpy (p->channel, hdr->channel_id, 3);
	p->channel[3] = '\0';
	strncpy (p->network, hdr->network_id, 2);
	strncpy (p->location, hdr->location_id, 2);
	p->network[2] = '\0';
 	p->sample_rate = hdr->sample_rate;
 	p->sample_rate_mult = hdr->sample_rate_mult;
	/* Assign desired gap and slew channel thresholds.		*/
	time_interval2 (1, p->sample_rate, p->sample_rate_mult, &seconds, &usecs);
	usecs += seconds * USECS_PER_SEC;
	p->gap_thresh = (gap_thresh >= 0) ? gap_thresh : 
	    (usecs / (USECS_PER_TICK*2));
	p->slew_thresh = (slew_thresh >= 0) ? slew_thresh : 
	    (usecs / (USECS_PER_TICK*2));

	/* Assign a type to the channel. */
	if (find_blockette(hdr,2000)) {
	    /* Opaque data channel. */
	    p->type = 'O';
	    p->ftype = 'O';
	}
	else if (p->sample_rate == 0 || p->sample_rate_mult == 0) {
	    /* Log channel. */
	    p->type = 'L';
	    p->ftype = 'L';
	}
	else {
	    /* Data channel. */
	    p->type = 'D';
	    p->ftype = 'D';
	}

	/* Keep a copy of the first header for the steam.		*/
	/* Zero any fields that will be used for channel totals.	*/
	p->init_hdr = dup_data_hdr(hdr);
	p->sum_hdr = dup_data_hdr(hdr);
	p->num_samples = hdr->num_samples;
	p->total_slew = 0;
	p->max_slew = 0;
	p->min_slew = 0;
	p->max_step = 0;
	p->num_blk = 0;
	p->prev_hdr = NULL;
	p->cur_hdr = hdr;
	p->data = p->databuffer;

	if ((! no_output) && xtractable(hdr)) {
	    sprintf (p->filename, "%s", p->station);
	    strcpy (p->filename, tempnam(".", p->filename));
	    if ( (p->fp = fopen(p->filename, "w+")) == NULL && errno == EMFILE) {
		/* Too many files.  Close a file unit and try again.	*/
		suspend_stream (&st_head);
		p->fp = fopen (p->filename, "w+");
	    }
	    if (p->fp == NULL) {
		fprintf (info, "unable to open file %s\n", p->filename);
		exit(1);
	    }
	    if (debug(DEBUG_STREAM)) {
		fprintf (info, "file %s open unit %d\n", p->filename, fileno(p->fp));
	    }
	    p->output = 1;
	}
	else {
	    p->output = 0;
	    p->fp = NULL;
	}
    }
    else {
	if (p->prev_hdr) free_data_hdr (p->prev_hdr);
	p->prev_hdr = p->cur_hdr;
	p->cur_hdr = NULL;
	p->cur_hdr = hdr;
	p->total_slew += dslew;
	if (p->total_slew > p->max_slew) p->max_slew = p->total_slew;
	if (p->total_slew < p->min_slew) p->min_slew = p->total_slew;
	if (fabs(dslew) > p->max_step) p->max_step = fabs(dslew);
	p->sum_hdr->num_samples += hdr->num_samples;
	p->sum_hdr->endtime = hdr->endtime;
    }

    /*	Check the sequence number to ensure that we have not skipped	*/
    /*	any blocks.							*/
    
    if (verbose && ignore_seq == 0 && seq_no != 0 && seq_no+1 != hdr->seq_no) {
	fprintf (info, "Warning: seq_no: expected %d, found %d\n",
		 seq_no+1, hdr->seq_no);
    }
    seq_no = hdr->seq_no;

    return (p);
}

/************************************************************************/
/*  write_channel:							*/
/*	Write component data to file.					*/
/************************************************************************/
int write_channel
   (ST_INFO	*st_p,		/* ptr to stream structure.		*/
    char	*buf,		/* ptr to sdr buffer.			*/
    int		blksize)	/* blksize.				*/
{
    int nd = blksize;
    int nw;
    char *data = buf;
    DATA_HDR *cur_hdr = st_p->cur_hdr;

    if (st_p->output == 0) return (0);
    if (st_p->fp == NULL && st_p->output) reopen_stream (st_p);

    /* Set sequence number in output seed record.	*/
    capnint(buf,st_p->num_blk,6);
    /* Add in network id field if necessary. */
    if (strlen(cur_hdr->network_id)==0 || strcmp(cur_hdr->network_id,"  ")==0) {
	capnstr(((SDR_HDR *)buf)->network_id,default_net,2);
    }
    while (nd > 0 && ( (nw = fwrite (data, sizeof(char), nd, st_p->fp)) != 0 ) ) {
	nd -= nw;
	data += nw;
    }
    if (nd != 0) {
	fprintf (info, "error writing file %s\n", st_p->filename);
	fprintf (info, "num_samples = %d\n", st_p->num_samples);
	exit(1);
    }

    return (st_p->cur_hdr->num_samples);
}

/************************************************************************/
/*  flush_channel:							*/
/*	Update the header block(s) for the channel(s).			*/
/************************************************************************/
int flush_channel
   (ST_INFO	*st_p,		/* ptr to stream structure.		*/
    int		flag)		/*  not used here.			*/
{
    char newname[1024], anewname[1024];
    SDR_TIME	btime;
    char	net[4];
    int		seconds, usecs;
    int		seq;
    INT_TIME	endtime;

    if (strlen(st_p->network)==0 || strcmp(st_p->network,"  ")==0) 
	strcpy(net,default_net); 
    else strcpy(net,st_p->network);
    trim(net);
    if (make_index) printf ("%s %s %s ", st_p->station, st_p->channel, net);
    if (gen_index) {
	time_interval2 (1, st_p->sum_hdr->sample_rate, st_p->sum_hdr->sample_rate_mult,
			&seconds, &usecs);
	endtime = add_time (st_p->sum_hdr->endtime, seconds, usecs);
	fprintf (stdout, "%s %s: rate=%d (%s to ",
		 st_p->sum_hdr->station_id, st_p->sum_hdr->channel_id, 
		 st_p->sum_hdr->sample_rate,
		 time_to_str(st_p->sum_hdr->begtime,date_fmt));
	fprintf (stdout, "%s) : %d points, %.1lf msec correction, ", 
		 time_to_str(endtime,date_fmt), st_p->sum_hdr->num_samples, 
		 (st_p->total_slew)/USECS_PER_MSEC);
	fprintf (stdout, "(min,max,max_step = %.1lf,%.1lf,%.1lf msec)\n",
		 (st_p->min_slew)/USECS_PER_MSEC, 
		 (st_p->max_slew)/USECS_PER_MSEC,
		 (st_p->max_step)/USECS_PER_MSEC);
    }
    if (st_p->output == 0) return (0);

    btime = encode_time_sdr (st_p->init_hdr->begtime,my_wordorder);
    if (passcal_filename) {
	if (split_on_rate) {
	    sprintf (newname, "%02d.%03d.%02d.%02d.%02d.%s.%s.%s.%s.%d",
		     btime.year%100, btime.day, btime.hour, btime.minute, 
		     btime.second, net, st_p->station, st_p->channel,
		     st_p->location, st_p->sample_rate);
	}
	else {
	    sprintf (newname, "%02d.%03d.%02d.%02d.%02d.%s.%s.%s.%s",
		     btime.year%100, btime.day, btime.hour, btime.minute, 
		     btime.second, net, st_p->station, st_p->channel,
		     st_p->location);
	}
    }
    else {
	if (split_on_rate) {
	    sprintf (newname, "%s.%s.%s.%s.%c.%d.%04d.%03d.%02d%02d%02d",
		     st_p->station, net, st_p->channel, st_p->location,
		     st_p->ftype,
		     st_p->sample_rate, btime.year, btime.day, 
		     btime.hour, btime.minute, btime.second);
	}
	else {
	    sprintf (newname, "%s.%s.%s.%s.%c.%04d.%03d.%02d%02d%02d",
		     st_p->station, net, st_p->channel, st_p->location,
		     st_p->ftype,
		     btime.year, btime.day, 
		     btime.hour, btime.minute, btime.second);
	}
    }
    seq = 0;
    strcpy (anewname, newname);
    if (unique_files) {
	while (access (anewname, F_OK) == 0) {
	    sprintf (anewname, "%s.%06d", newname, seq++);
	    if (seq >= 1000000) {
		fprintf (info, "Error: ran out of sequence numbers for %s\n", newname);
		exit(1);
	    }
	}
    }
	
    if (rename (st_p->filename, anewname) != 0) {
	fprintf (stderr, "unable to rename %s to %s\n",
		 st_p->filename, anewname);
	return (EOF);
    }
    else strcpy (st_p->filename, newname);
    if (st_p->fp != NULL) {
	if (debug(DEBUG_STREAM)) {
	    fprintf (info, "file %s close unit %d\n", st_p->filename, fileno(st_p->fp));
	}
	fclose (st_p->fp);
    }
    return (0);
}

/************************************************************************/
/*  decode_hdr:								*/
/*	Get information about this block of data.			*/
/************************************************************************/
DATA_HDR *decode_hdr
   (char	*raw_hdr,	/* ptr to raw block header.		*/
    int		maxbytes)	/* max number of bytes in block.	*/
{
    DATA_HDR	*hdr;
    int		type;

    type = get_header_type (raw_hdr, maxbytes);
    switch (type) {
    case SDR_HDR_TYPE:
	hdr = decode_hdr_sdr ((SDR_HDR *)raw_hdr, maxbytes);
	break;
    case SDR_VOL_HDR_TYPE:
	hdr = decode_hdr_sdr ((SDR_HDR *)raw_hdr, maxbytes);
	break;
    default:
	fprintf (info, "unknown header type: %d\n", type);
	exit(1);
	break;
    }
    if ( (hdr == NULL) && (herrno != 0) ) {
	fprintf (info, "Invalid header found\n");
	exit(1);
    }

    return (hdr);
}

/************************************************************************/
/*  get_header_type:							*/
/*	Determine the type of input stream.				*/
/************************************************************************/
int get_header_type 
   (char    *p,			/* ptr to buffer containing header.	*/
    int	    nr)			/* max number of bytes for header.	*/
{
    if (is_drm_header((STORE_FILE_HEAD *)p,nr)) return DRM_HDR_TYPE;
    else if (is_qda_header((QDA_HDR *)p,nr)) return QDA_HDR_TYPE;
    else if (is_sdr_header((SDR_HDR *)p,nr)) return SDR_HDR_TYPE;
    else if (is_sdr_vol_header((SDR_HDR *)p,nr)) return (SDR_VOL_HDR_TYPE);
    else return (-1);
}

/************************************************************************/
/*  suspend_stream:							*/
/*	Find an open stream, save its position, and close it.		*/
/*  return:								*/
/*	0 on success, -1 on failure.					*/
/************************************************************************/
int suspend_stream 
   (ST_INFO	*phead)		/* ptr to stream struct linked list.	*/
{
    static int warning_sent = 0;
    ST_INFO *st_p;
    for (st_p = phead->next; st_p != NULL; st_p = st_p->next) {
	if (st_p->fp != NULL && st_p->fp != stderr &&
	    st_p->fp != stdin && st_p->fp != stdout &&
	    (st_p->pos = ftell (st_p->fp)) >= 0) {
	    if (debug(DEBUG_STREAM) && ! warning_sent) {
		fprintf (info, "Warning: too many open files -- reclaiming file units\n"),
		warning_sent = 1;
	    }
	    if (debug(DEBUG_STREAM)) {
		fprintf (info, "suspend stream: %s\n", st_p->filename);
	    }
	    fclose (st_p->fp);
	    st_p->fp = NULL;
	    return (0);
	}
    }
    return (-1);
}

/************************************************************************/
/*  reopen_stream:							*/
/*	Reopen stream that was temporarily closed.			*/
/************************************************************************/
int reopen_stream
   (ST_INFO	*st_p)		/* ptr to stream struct to reopen.	*/
{
    st_p->fp = fopen (st_p->filename, "a");
    if (st_p->fp == NULL) {
	suspend_stream (&st_head);
	st_p->fp = fopen (st_p->filename, "a");
    }
    if (st_p->fp == NULL) {
	fprintf (info, "unable to re-open file %s\n", st_p->filename);
	exit (1);
    }
    if (debug(DEBUG_STREAM)) {
	fprintf (info, "reopen_stream:  %s unit %d\n", st_p->filename, fileno(st_p->fp));
    }
    /*::
    if (fseek (st_p->fp, st_p->pos, SEEK_SET) != 0) {
	fprintf (info, "unable to re-position file %s\n", st_p->filename);
	exit (1);
    }
    ::*/
    return (0);
}

/************************************************************************/
/*  parse_chanlist:							*/
/*	Parse channel extraction list.					*/
/************************************************************************/
int parse_chanlist (char *list)
{
    char entry[20];
    char *p, *p0, *p1;
    int l;

    xtractc = 0;
    xtractv = (CHANLIST **)malloc(sizeof(CHANLIST *));
    if (xtractv == NULL) {
	fprintf (stderr, "Error mallocing xtractv\n");
	exit (1);
    }
    xtractv[xtractc] = NULL;
    p = list;
    while (p != NULL && *p != '\0') {
	p1 = strchr (p,',');
	if (p1 == NULL) p1 = p+strlen(p);
	l = p1-p;
	strncpy (entry, p, l);
	entry[l] = '\0';
	p = (*p1) ? p1+1 : p1;
	xtractv[xtractc] = (CHANLIST *)malloc(sizeof(CHANLIST));
	memset((void *)xtractv[xtractc], 0, sizeof(CHANLIST));
	strcpy(xtractv[xtractc]->channel,"*");
	strcpy(xtractv[xtractc]->network,"*");
	strcpy(xtractv[xtractc]->location,"*");

	/* A single componet:	station					*/
	/* 2 components:	station.channel				*/
	/* 3 components:	station.network.channel			*/
	/* 4 components:	station.network.channel.location.	*/
	p0 = &entry[0];
	p1 = strchr(p0,'.');
	if (p1 == NULL) {
	    strcpy(xtractv[xtractc]->station,p0);
	}
	else {
	    /* station is always first component.			*/
	    l=p1-p0;
	    strncpy(xtractv[xtractc]->station,p0,l);
	    xtractv[xtractc]->station[l] = '\0';
	    p0 = p1+1;
	    /* 2 components is station.channel.				*/
	    /* More than 2 components is station.net.channel...		*/
	    p1 = strchr(p0,'.');
	    if (p1 != NULL) {
		l=p1-p0;
		strncpy(xtractv[xtractc]->network,p0,l);
		xtractv[xtractc]->network[l] = '\0';
		p0 = p1+1;
	    }
	    /* Next component is channel.				*/
	    p1 = strchr(p0,'.');
	    if (p1 != NULL) {
		l=p1-p0;
		strncpy(xtractv[xtractc]->channel,p0,l);
		xtractv[xtractc]->channel[l] = '\0';
		strcpy(xtractv[xtractc]->location,p1);
	    }
	    else {
		strcpy(xtractv[xtractc]->channel,p0);
	    }
	}
	uppercase(xtractv[xtractc]->station);
	uppercase(xtractv[xtractc]->network);
	uppercase(xtractv[xtractc]->channel);
	uppercase(xtractv[xtractc]->location);
	xtractc++;
	xtractv = (CHANLIST **)realloc(xtractv, (xtractc+1)*sizeof(CHANLIST *));
	if (xtractv == NULL) {
	    fprintf (stderr, "Error reallocing xtractv\n");
	    exit (1);
	}
	xtractv[xtractc] = NULL;
    }
    return (xtractc);
}

int xtractable (DATA_HDR *hdr) 
{
    int i;
    if (xtractc == 0) return (1);

    for (i=0; i<xtractc; i++) {
	if (fnmatch (xtractv[i]->station, hdr->station_id, 0) == 0 &&
	    fnmatch (xtractv[i]->network, hdr->network_id, 0) == 0 &&
	    fnmatch (xtractv[i]->channel, hdr->channel_id, 0) == 0 &&
	    fnmatch (xtractv[i]->location, hdr->location_id, 0) == 0)
	    return (1);
    }
    return (0);
}

int different_day (INT_TIME begtime1, INT_TIME begtime2)
{
    EXT_TIME et1;
    EXT_TIME et2;
    et1 = int_to_ext(begtime1);
    et2 = int_to_ext(begtime2);
    return (et1.year != et2.year || et1.doy != et2.doy);
}

int first_b2000_seq (DATA_HDR *hdr)
{
    BS *bs = hdr->pblockettes;
    int seq;
    bs = find_pblockette (hdr, bs, 2000);
    if (bs == NULL) return (-1);
    seq = ((BLOCKETTE_2000 *)(bs->pb))->record_num;
    return (seq);
}

int last_b2000_seq (DATA_HDR *hdr)
{
    BS *bs = hdr->pblockettes;
    BS *nextbs = NULL;
    int seq;
    bs = find_pblockette (hdr, bs, 2000);
    if (bs == NULL) return (-1);
    nextbs = bs->next;
    while ((nextbs = find_pblockette (hdr, nextbs, 2000)) != NULL) {
	bs = nextbs;
	nextbs = bs->next;
    }
    seq = ((BLOCKETTE_2000 *)(bs->pb))->record_num;
    return (seq);
}
