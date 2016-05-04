/************************************************************************/
/*  Stream information structure.  One allocated for each stream.	*/
/*  This is used to store stream information in a manner that is	*/
/*  independent of both the input and output format.			*/
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

/*	$Id: st_info.h,v 1.8 2013/09/17 23:47:02 doug Exp $ 	*/

#include    "qlib2.h"

typedef struct _st_info {
    char	    station[6]; /* name of station.			*/
    char	    channel[4];	/* name of component.			*/
    char	    network[3];	/* name of network.			*/
    char	    location[3];/* name of location.			*/
    char	    type;	/* type of file.			*/
    char	    ftype;	/* type of label for file.		*/
    char	    dir[256];	/* directory for files.			*/
    char	    filename[256];  /* filename of output data file.	*/
    char	    tfile[256];	/* filename of time corr data file.	*/
    int		    sample_rate; /* data rate.				*/
    int		    sample_rate_mult; /* data rate_mult.		*/
    int		    num_blk;	/* # of data block read for this comp.	*/
    double	    cur_slew;	/* current slew for this block.		*/
    double	    total_slew;	/* total time slew in ticks.		*/
    double	    max_slew;	/* maximum slew at any time.		*/
    double	    min_slew;	/* minimum slew at any time.		*/
    double	    max_step;	/* abs(max slew step per block).	*/
    int		    gap_thresh;	/* threshold for inter-block gap.	*/
    int		    slew_thresh;/* threshold for total channel slew.	*/
    int		    output;	/* boolean flag to output data.		*/
    long int	    pos;	/* position in output file.		*/
    FILE	    *fp;	/* file pointer for output data file.	*/
    FILE	    *tp;	/* file pointer for time corr data file.*/
    DATA_HDR	    *init_hdr;	/* ptr to initial header structure.	*/
    DATA_HDR	    *sum_hdr;	/* ptr to summary header structure.	*/
    DATA_HDR	    *prev_hdr;	/* ptr to previous header structure.	*/
    DATA_HDR	    *cur_hdr;	/* ptr to current header structure.	*/
    int		    num_samples;/* # of samples, pointed to by data.	*/
    int		    *data;	/* ptr to start of decompressed data.	*/
    int		    *databuffer;/* buffer for decompressed data.	*/
    struct _st_info *next;	/* ptr to next stream structure in list.*/
    struct _st_info *prev;	/* ptr to prev stream structure in list.*/
} ST_INFO;

