/************************************************************************/
/*  Routines for sdrsplit.						*/
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

/*	$Id: procs.h,v 1.4 2010/03/03 03:53:01 doug Exp $ 	*/

int blkread 
   (FILE	*in,		/* input FILE.				*/
    char	*block,		/* sdr record.				*/
    int		blksize);	/* # of bytes to read (0 -> MAX_BLKSIZE)*/

ST_INFO *alloc_channel (void);

void free_channel
   (ST_INFO	*p);		/* ptr to stream structure.		*/

ST_INFO *find_st
   (DATA_HDR	*hdr);		/* ptr to DATA_HDR for record.		*/

int write_channel
   (ST_INFO	*st_p,		/* ptr to stream structure.		*/
    char	*buf,		/* ptr to sdr buffer.			*/
    int		blksize);	/* blksize.				*/

int flush_channel
   (ST_INFO	*st_p,		/* ptr to stream structure.		*/
    int		flag);		/*  not used here.			*/

DATA_HDR *decode_hdr
   (char	*raw_hdr,	/* ptr to raw block header.		*/
    int		maxbytes);	/* max number of bytes in block.	*/

int get_header_type 
   (char    *p,			/* ptr to buffer containing header.	*/
    int	    nr);		/* max number of bytes for header.	*/

int parse_chanlist (char *list);
