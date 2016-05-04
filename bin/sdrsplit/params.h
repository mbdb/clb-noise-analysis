/************************************************************************/
/*  Parameters for sdrsplit.						*/
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

/*	$Id: params.h,v 1.4 2000/10/26 13:03:47 doug Exp $ 	*/

#define	DEBUG
#define	DEBUG_SLEW	1
#define	DEBUG_TIME	2
#define	DEBUG_BLOCK	4
#define	DEBUG_HDRS	8
#define	DEBUG_FLAG	16
#define	DEBUG_STREAM	32
#define	DEBUG_OUT	64

#define	DEBUG_RATE	1024	/* Do not include in DEBUG_ANY	*/

#define	DEBUG_ANY	\
    (DEBUG_BLOCK | DEBUG_HDRS | DEBUG_SLEW | DEBUG_TIME | DEBUG_FLAG | \
     DEBUG_STREAM | DEBUG_OUT)
#define	debug(val)  (debug_option & (val))

/*  Keys used to control stream header creation and manipulation.	*/
#define	HDR_INIT	1
#define	HDR_FLUSH	2

#define	MAX_BLKSIZE	8192

/* Define relation operators for blockette 200 record numbers.		*/
/* HR2000 is half of the range for blockette 2000 record numbers.	*/
#define	MAX2000	0x7FFFFFFF
#define HR2000	0x40000000 
#define	EQ2000(a,b) ( a==b )
#define	GT2000(a,b) ( ((a>b) && (a-b)<=HR2000) || ((a<b) && (b-a) > HR2000) )
#define	GE2000(a,b) ( EQ2000(a,b) || GT2000(a,b) )
#define	NE2000(a,b) ( ! EQ2000(a,b) )
#define	LT2000(a,b) ( GT2000(b,a) )
#define	LE2000(a,b) ( GE2000(b,a) )
#define	NEXT2000(a) ( (a<0) ? a : (a>=MAX2000) ? 0 : a+1 )
#define	PREV2000(a) ( (a<0) ? a : (a==0) ? MAX2000 : a-1 )
#define	DIFF2000(a,b) (EQ2000(a,b) ? 0 : \
		       GT2000(a,b) ? (((a-b)+HR2000)%HR2000) : \
		       (-1*(((b-a)+HR2000)%HR2000)))
#define IGNORE2000  -1
#define	DEFAULT_BLOCK_TOL   100
