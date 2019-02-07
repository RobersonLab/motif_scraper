import argparse
import logging
import multiprocessing as mp
import os
import pyfaidx
import regex
import six
import shutil
import sys
import time

from motif_scraper import __script_name__, __version__, rev_comp, make_degenerate_regex, SequenceMotif, fasta_motif_scan

def main():
	#############
	# arg parse #
	#############
	parser = argparse.ArgumentParser( prog=__script_name__, epilog="%s v%s" % ( __script_name__, __version__ ) )

	# FASTA index will be created if it does not exist when pyfaidx Fasta is initialized
	parser.add_argument( 'fastaFile', help="Path to FASTA file to be scanned for degenerate sequence motifs." )
	parser.add_argument( '--motif', '-m', help="A degenerate sequence motif. Can be specified multiple times.", action='append' )
	parser.add_argument( '--motif_file', help="A file containing motifs to search file. Can be specified multiple times.", action='append' )
	parser.add_argument( '--outputFile', '-o', help="Defaults to detected_motifs.csv", default="detected_motifs.csv" )
	parser.add_argument( '--region', '-r', help="Defines region to search for sites. Use 'contig:start-end' for regions, or 'contig' for whole contig. If no regions are specified, the entire FASTA file will be searched! Starts expected to be *0 start* start and Ends *1 start*", action='append', required=False )
	parser.add_argument( '--cores', '-p', help="Run search on multiple contigs / strands simultaneously", type=int, default=1 )
	parser.add_argument( '--valid_regex', help="Query is valid regex. *WILL NOT* reverse complement. Specify sequence and strand with careful consideration.", default=False, action='store_true' )
	# molecule
	parser.add_argument( '--search_strand', '-s', help="Default searches both strands, but can be set to only one.", choices=[ '+1', '-1', 'both' ], default='both' )
	parser.add_argument( '--no_motif_overlaps', help="If this option is set it turns off overlapping motif matches.", default=True, action='store_false' )
	parser.add_argument( '--file_buffer', help="On low memory systems a file buffer can be used, which is slower by has low memory requirements.", default=False, action='store_true' )
	parser.add_argument( '--copy_buffer_size', help="Size of buffer (MB) for file copy", type=int, default=10 )
	parser.add_argument( "--loglevel", choices=[ 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL' ], default='INFO' )

	args = parser.parse_args()
	
	#################
	# setup logging #
	#################
	logging.basicConfig( format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s' )
	logger = logging.getLogger( __script_name__ )
	logger.setLevel( args.loglevel )
	
	########################
	# fix up input parsing #
	########################
	allow_overlapping_motifs = not args.no_motif_overlaps
	use_file_buffer = args.file_buffer
	file_buffer_size = args.copy_buffer_size * 1024 * 1024
	
	if args.search_strand == "both":
		run_strands = ( '+', '-' )
	elif args.search_strand == '-1':
		run_strands = ( '-', )
	else:
		run_strands = ( '+', )
	
	if args.motif is None and args.motif_file is None:
		logger.critical( "Must specify at least one --motif or path to a --motif_file" )
		sys.exit( 1 )
		
	if args.motif is None:
		args.motif = []
	
	if args.motif_file is not None:
		for path in args.motif_file:
			with open( path, 'r' ) as MOTIF_IN:
				for line in MOTIF_IN:
					line = line.strip()
					
					if len( line ) == 0:
						continue
					elif line[0] == '#':
						continue
					args.motif.append( line.upper() )
				
	args.motif = list( set( [ seq.upper() for seq in args.motif ] ) ) # uniquify motifs in case they show up multiple times
	
	if len( args.motif ) == 0:
		logger.critical( "No motifs specified!" )
		sys.exit( 1 )
	
	######################
	# set options string #
	######################
	options = "\n%s v%s\n\nOptions\n=======\n" % ( __script_name__, __version__ )
	options += "FASTA: %s\n" % ( args.fastaFile )
	options += "Motifs to search: %s\n" % ( str( args.motif ) )
	options += "Strands of FASTA to search: %s\n" % ( str( run_strands ) )
	options += "Output file: %s\n" % ( args.outputFile )
	options += "Regions: %s\n" % ( "All contigs" if args.region is None else str( args.region ) )
	options += "Motifs are valid regex (no processing): %s\n" % ( str( args.valid_regex ) )
	options += "Allow overlapping motifs?: %s\n" % ( str( allow_overlapping_motifs ) )
	options += "Processes: %s\n" % ( args.cores ) # cores
	options += "Buffering: %s\n" % ( "Memory" if not args.file_buffer else "File" )
	if args.file_buffer:
		options += "File copy buffer size (MB): %s\n" % ( args.copy_buffer_size )
	options += "Log level: %s\n" % ( str( args.loglevel ) )
	
	logger.info( options )
		
	###################
	# Open FASTA file #
	###################
	region_list = []
	
	with pyfaidx.Fasta( args.fastaFile, as_raw=True ) as FAIDX:
		if args.region is None:
			args.region = FAIDX.keys()
		
		for reg in args.region:
			contig, start, end = pyfaidx.ucsc_split( reg )  # returns [0,1) coordinates with NoneType if not start or end

			# no need for checking for contig in FASTA file, as this is done by pyfaidx
			if start is None:
				start = 0 # pyfaidx using 0-base start. Remember for custom regions!!!
			if end is None:
				end = len( FAIDX[contig] )
			
			for curr_motif in args.motif:
				for curr_strand in run_strands:
					region_list.append( ( curr_motif, contig, start, end, curr_strand ) )
					
					logger.debug( "%s\n" % ( str( region_list[-1] ) ) )

	################
	# process data #
	#############################
	# async process the contigs #
	#############################
	work_pool = mp.Pool( processes=args.cores )
	working_data = [ work_pool.apply_async( fasta_motif_scan, args=( args.fastaFile, x, args.valid_regex, allow_overlapping_motifs, use_file_buffer ) ) for x in region_list ]
	work_output = [ x.get() for x in working_data ]
	
	##############################
	# figure order out           #
	# b/c executed asyncronously #
	##############################
	site_count = 0
	
	result_order = []
	
	for reg in region_list:
		index = 0
		
		for out, dict, fname, current_site_count in work_output:
			if out == reg:
				logger.debug( "Returned result %s matches %s" % ( index, ' '.join( [ str(x) for x in reg ] ) ) )
				
				result_order.append( index )
				site_count += current_site_count
				
				break
			else:
				index += 1
				
	if len( result_order ) != len( region_list ):
		logger.critical( "Did not find results for all requested regions" )
		sys.exit( 1 )
	
	#################
	# write outputs #
	#################
	if use_file_buffer == True:
		with open( args.outputFile, 'wb' ) as OUTFH:
			# self.contig, self.start, self.end, self.strand, self.seq, self.motif
			OUTFH.write( "Contig,Start,End,Strand,Sequence,Motif\n" )
			
			for idx in result_order:
				curr_fname = work_output[idx][2]
				with open( curr_fname, 'rb' ) as COPYINPUT:
					shutil.copyfileobj( COPYINPUT, OUTFH, file_buffer_size )
				os.remove( curr_fname )
	else:
		with open( args.outputFile, 'w' ) as OUTFH:
			# self.contig, self.start, self.end, self.strand, self.seq, self.motif
			OUTFH.write( "Contig,Start,End,Strand,Sequence,Motif\n" )
			
			for idx in result_order:
				for site in work_output[idx][1]:
					OUTFH.write( "%s\n" % ( site ) )

	logger.info( "%s total sites found" % ( site_count ) )

if __name__ == "__main__":
	main()
	