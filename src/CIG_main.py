"""
Copyright (c) <2017> <Dr. Kaifu Chen lab, Research Institute, Houston Methodist Hospital >

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import argparse, sys, os
from queryUtils import GEO_query
from search import SOFTQuickParser
from Related_Sample_Search import Related_Sample_Search
from setup import get_settings, setup
from update import update


def Help():
    print "\nDanpos Grid Optimization"
    print "A list of functions for grid optimization, please try:\npython danpos.py grid -h"
    print "\nFuctions:"
    print "\tgrid:\n\tfind the best parameter for Danpos."
    print ""

def CIG_grid():
    '''
    this function provie an entrance to search function

    '''
    if (len(sys.argv) < 3) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least one parameter need to be specified, will print help message if no parameter is specified
        print "\nusage:\n\npython danpos.py grid [optional arguments] <target_table1> <target_table2> <danpos_result_path> <features> <output_prefix> <output_path>\n\nfor more help, please try: python danpos.py grid -h\n"
        return 1

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\npython danpos.py grid [optional arguments] <target_table1> <target_table2> <danpos_result_path> <features> <output_prefix> <output_path>\n\n",
                                     description='',epilog="Chen lab, Houston Methodist")
    parser.add_argument('command', default=None, help="set as 'grid' to perform parameter optimization")

    parser.add_argument('target_table1', default=None,
                        help="The first table of genes, containing the columns at least 'gene_id', 'sample_prefix', table delimiter is recognized by the file surfix (csv, tsv, txt, xls, xlsx)", )
    parser.add_argument('target_table2', default=None,
                        help="The second table of genes, containing the columns at least 'gene_id', 'sample_prefix', table delimiter is recognized by the file surfix (csv, tsv, txt, xls, xlsx)")
    parser.add_argument('danpos_result_path', default=None,
                        help="folder containing the danpos peak calling result tables, make sure the tables startswith sample_prefix and with '_' as delimiter in the name")
    parser.add_argument('features', default=None,
                        help="")
    parser.add_argument('output_prefix', default=None,
                        help="the prefix for output files")
    parser.add_argument('output_path', default=None,
                        help="the output path")
    parser.add_argument('up_stream_grid', default=None,
                        help="the optimization grid for upstream distance, in the format:'10000:2:1000', meaning, the start grid is 10000, every iteration the grid shrink for 2 times, and the final grid need to be larger than 1000")
    parser.add_argument('down_stream_grid', default=None,
                        help="the optimization grid for downstream distance, in the format:'10000:2:1000', meaning, the start grid is 10000, every iteration the grid shrink for 2 times, and the final grid need to be larger than 1000")
    parser.add_argument('height_grid', default=None,
                        help="the optimization grid for height, in the format:'10000:2:1000', meaning, the start grid is 10000, every iteration the grid shrink for 2 times, and the final grid need to be larger than 1000")

    ## optional parameters
    parser.add_argument('-f', dest='function', metavar='', default='wilcoxon',
                        help="wilcoxon or fisher")

    args = None

    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print "\Chipseqpair\n"
        parser.print_help()
        print "\n"
        return 0
    elif len(sys.argv) >= 3:
        try:
            args = parser.parse_args()
        except:
            print "\nfor more help, please try: python CSP.py search -h\n"
            return 1

    if args is not None:
        settings = get_settings()
        encode_pkl = settings['Encode']
        roadmap_pkl = settings['Roadmap']
        GGRmap_pkl = settings['GGR']
        GSMGSE_pkl = settings['GSMGSE_pkl_path']

        keywords = args.feature_key_words.split(",")

        output_prefix = args.output_prefix
        if output_prefix is None:
            output_prefix = keywords[0]

        output_path = args.output_path
        if output_path is None:
            output_path = './search_output/'
            if not os.path.isdir(output_path):
                os.system("mkdir search_output")

        if args.keywords_begin == '':
            keywords_begin = []
        else:
            keywords_begin = args.keywords_begin.split(",")

        type_seq = args.type_seq
        ignorcase = args.ignorecase
        geo = args.geo
        geo_file = args.geo_file

        species = args.species
        encode_remove = 1
        roadmap_remove = 1

        cwd = args.MetaData
        process = args.process

        if cwd is None:
            cwd = settings['MetaData']

        if cwd is None or cwd == "None":
            cwd = None
            encode_remove = True
            roadmap_remove = True

        SOFTQuickParser(output_prefix, output_path, keywords, keywords_begin,
                        type_seq=type_seq, ignorecase=ignorcase, geo=geo, geofile=geo_file, output_type=species,
                        encode_remove=encode_remove, roadmap_remove=roadmap_remove,
                        encode_pkl=encode_pkl, roadmap_pkl=roadmap_pkl, GGRmap_pkl=GGRmap_pkl,
                        GSMGSE_pkl=GSMGSE_pkl, cwd=cwd, process=process)
        return
    return 1

if len(sys.argv) > 1:
    if sys.argv[1] == "grid":
        CIG_grid()
    else:
        Help()
else:
    print "\nGrid Optimization"
    print "A list of functions for grid optimization, please try:\npython grid.py -h"
    print "\nFuctions:"
    print "\tgrid:\n\tfind the best parameter for chipseq parameters."
    print ""
