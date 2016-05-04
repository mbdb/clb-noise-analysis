/************************************************************************/
/*  Version number.							*/
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

/*	$Id: version.h,v 1.13 2013/09/17 23:47:02 doug Exp $ 	*/

#define	VERSION	"1.4.12 (2013.260)"

/*
Modification History
------------------------------------------------------------------------
Date	    Ver			Who	What
------------------------------------------------------------------------
2013/09/17  1.4.12 (2013.260)	Doug Neuhauser
	Fixed most compiler warning msgs, updated Makefile.
2010/07/16  1.4.11		Doug Neuhauser
	Updated to merge contiguous segments of data per channel that are 
	separated in input file by other non-contiguous data of the same 
	channel.
2010/04/26  1.4.10		Doug Neuhauser
	Link with updated qlib2 to support 'M' MiniSEED records.
2010/03/16  1.4.9		Doug Neuhauser
	Fixed eof detection (again).
2010/03/09  1.4.8		Doug Neuhauser
	Fixed eof detection.
2010/03/01  1.4.7		Doug Neuhauser
	Updated to with with multiple blksizes in input file.
2008/02/22  1.4.6		Doug Neuhauser
	Added -u option to create unique files.
2004/12/11  1.4.5		Doug Neuhauser
	Use network in comparison.
2003/01/05  1.4.4
	Use qlib2 function is_data_hdr_ind() to accept either R,D,Q in record.
2000/10/26  1.4.3 (2000.300)	Doug Neuhauser
	sdrsplit.c, params.h,st_info.h, version.h
	Program supplied appropriate filetype character (D,L,O) 
	based on contents of record.
1999/04/07  1.4.2 (1999.098)	Doug Neuhauser
	sdrsplit.c
	fnmatch.c
	Converted to POSIX and XPG4 fnmatch, which returns 0 on success.
*/
