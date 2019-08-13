#! /usr/bin/python
#########################################################################################
#                                                                                       #
# This program is free software: you can redistribute it and/or modify it under the     #
# terms of the GNU General Public License as published by the Free Software Foundation, #
# either version 3 of the License, or (at your option) any later version.               #
#                                                                                       #
#########################################################################################
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, ##
## INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A       ##
## PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT  ##
## HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION   ##
## OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE      ##
## SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                              ##
#########################################################################################

import os
import sys

import logging
import subprocess
import re 

from biolib.common import (make_sure_path_exists, 
                            check_dir_exists,
                            check_file_exists,
                            remove_extension)

from mama.common import findEx
from mama.defaultValues import DefaultValues


class OptionsParser():

    def __init__(self):
        """Initialization"""
        self.logger = logging.getLogger('timestamp')

    def features(self,options):
        """Making bam features matrix"""

        make_sure_path_exists(options.output_dir)
        reads_abundance = os.path.join(options.output_dir,  DefaultValues.FEATURES_ABUNDANCE_FILES[0])
        reads_normalised = os.path.join(options.output_dir, DefaultValues.FEATURES_ABUNDANCE_FILES[1])
        reads_relative = os.path.join(options.output_dir,   DefaultValues.FEATURES_ABUNDANCE_FILES[2])
        base_abundance = os.path.join(options.output_dir,   DefaultValues.FEATURES_ABUNDANCE_FILES[3])
        base_normalised = os.path.join(options.output_dir,  DefaultValues.FEATURES_ABUNDANCE_FILES[4])
        base_relative = os.path.join(options.output_dir,    DefaultValues.FEATURES_ABUNDANCE_FILES[5])

        features_size = {}
        counts = {}
        counts_base = {}

        self.logger.info('Get features and initialise matrix')
        with open(options.faidx) as f:
            for line in f:
                if not line.startswith('#') :
                    line_list = line.rstrip().split('\t')
                    features = line_list[0]
                    features_size[features] = line_list[1]
                    counts[features] = 0 
                    counts_base[features] = 0 

        counts_all = []
        counts_all_normalised = []
        counts_all_relative = []
        counts_base_all = []
        counts_base_all_normalised = []
        counts_base_all_relative = []

        header = ["Features","Features_size"]
        self.logger.info('Browse alignement file(s)')

        samtoolsexec=findEx('samtools')
        samtoolsthreads='-@ ' + options.threads
        samtoolsminqual='-q ' + options.mapQ

        with open(options.bam_list,'r') as b :
            for bam in b :
                if bam.startswith('#') :
                    continue
                i = 0
                alignementfile,librarysize = bam.split('\t')
                if librarysize == '' or librarysize == 0 :
                    librarysize = 1 
                samplename = remove_extension(os.path.basename(alignementfile),options.extension)
                header.append(samplename)
                self.logger.info('\t'+samplename)
                cmd = [ samtoolsexec,'view',samtoolsthreads,samtoolsminqual, alignementfile ]
                p = subprocess.Popen(cmd, stdout=subprocess.PIPE).stdout
                for line in p:
                    line=line.decode(sys.getdefaultencoding()).rstrip()
                    if i > 0 and i % 1000000 == 0 :
                        self.logger.info("Alignment record %s processed" % i )
                    i += 1
                    line_list = line.split('\t')
                    features = line_list[2]
                    cigar = line_list[5]
                    base_mapped = 0 
                    match = re.findall(r'(\d+)M',cigar)
                    read_len = len(line_list[6])  
                    for base_match in match :
                        base_mapped += int(base_match)
                    if read_len == 0 :
                        self.logger.info(line_list)

                    if base_mapped/read_len < float(options.id_cutoff) : 
                        continue 
                
                    counts[features] += 1

                    if options.discard_gene_length_normalisation :
                        counts_base[features] += base_mapped
                    else :
                        counts_base[features] += base_mapped / int(features_size[features])
                
                # raw reads count 
                counts_all.append(counts.copy())
            
                # normalised reads count
                count_tmp = {}
                count_tmp = {k: (v / int(librarysize))*options.feature_normalisation for k, v in counts.items()} 
                counts_all_normalised.append(count_tmp.copy())

                # relative reads count
                count_tmp = {}
                count_tmp = {k: v / total for total in (sum(counts.values()),) for k, v in counts.items()}
                counts_all_relative.append(count_tmp.copy())

                # raw bases count 
                counts_base_all.append(counts_base.copy())

                # normalised bases count
                count_tmp = {}
                count_tmp = {k: (v / int(librarysize))*options.feature_normalisation for k, v in counts_base.items()} 
                counts_base_all_normalised.append(count_tmp.copy())

                # relative bases count
                count_tmp = {}
                count_tmp = {k: v / total for total in (sum(counts_base.values()),) for k, v in counts_base.items()} 
                counts_base_all_relative.append(count_tmp.copy())

                for fn in counts:
                    counts[fn] = 0
                    counts_base[fn] = 0

        self.logger.info('Print matrices')
    
        self.logger.info('Print raw reads abundance matrix in %s' % reads_abundance)
        output_handle = open(reads_abundance, "w")
        output_handle.write('\t'.join(header)+'\n')
        for fn in counts.keys():
            if sum([c[fn] for c in counts_all]) == 0 and options.removed :
                continue
            else :
                output_handle.write('\t'.join([fn] + [features_size[fn]] + [str(c[fn]) for c in counts_all])+'\n') 
        output_handle.close()
        
        self.logger.info('Print normalised reads abundance matrix in %s' % reads_normalised)
        output_handle = open(reads_normalised, "w")
        output_handle.write('\t'.join(header)+'\n')
        for fn in counts.keys():
            if sum([c[fn] for c in counts_all_normalised]) == 0 and options.removed :
                continue
            else :
                output_handle.write('\t'.join([fn] + [features_size[fn]] + [str(c[fn]) for c in counts_all_normalised])+'\n') 
        output_handle.close()

        self.logger.info('Print relative reads abundance matrix in %s' % reads_normalised)
        output_handle = open(reads_relative, "w")
        output_handle.write('\t'.join(header)+'\n')
        for fn in counts.keys():
            if sum([c[fn] for c in counts_all_relative]) == 0 and options.removed :
                continue
            else :
                output_handle.write('\t'.join([fn] + [features_size[fn]] + [str(c[fn]) for c in counts_all_relative])+'\n')
        output_handle.close()

        self.logger.info('Print raw base abundance matrix in %s' % reads_normalised)
        output_handle = open(base_abundance, "w")
        output_handle.write('\t'.join(header)+'\n')
        for fn in counts_base.keys():
            if sum([c[fn] for c in counts_all]) == 0 and options.removed :
                continue
            else :
                output_handle.write('\t'.join([fn] + [features_size[fn]] + [str(c[fn]) for c in counts_base_all])+'\n') 
        output_handle.close()

        self.logger.info('Print normalised base abundance matrix in %s' % reads_normalised)
        output_handle = open(base_normalised, "w")
        output_handle.write('\t'.join(header)+'\n')
        for fn in counts_base.keys():
            if sum([c[fn] for c in counts_all_normalised]) == 0 and options.removed :
                continue
            else :
                output_handle.write('\t'.join([fn] + [features_size[fn]] + [str(c[fn]) for c in counts_base_all_normalised])+'\n') 
        output_handle.close()

        self.logger.info('Print relative base abundance matrix in %s' % reads_normalised)
        output_handle = open(base_relative, "w")
        output_handle.write('\t'.join(header)+'\n')
        for fn in counts_base.keys():
            if sum([c[fn] for c in counts_all_relative]) == 0 and options.removed :
                continue
            else :
                output_handle.write('\t'.join([fn] + [features_size[fn]] + [str(c[fn]) for c in counts_base_all_relative])+'\n')
        output_handle.close()

        self.logger.info('Matrices printed')


    def annoted_features(self,options):
        """Making annoted features matrix"""

        missing = []

        features2annotation = {}
        with open(options.features_annotation) as f:
            for line in f: 
                line = line.rstrip()
                features_id,annotation = line.split('\t')
                features2annotation[features_id] = annotation
        
        counts = {}
        id2description = {}
        annotation_id_list = []
        with open(options.annotation_description) as f:
            for line in f:
                line = line.rstrip()
                annotation_id,description = line.split('\t')
                id2description[annotation_id] = description                
                annotation_id_list.append(annotation_id) 
                counts[annotation_id] = {} 

        
        annotation_id_list.append('hypothetical protein')
        counts['hypothetical protein'] = {} 

        check_dir_exists(options.features_dir)        
        input_matrices = DefaultValues.FEATURES_ABUNDANCE_FILES
        output_matrices = DefaultValues.ANNOTATE_ABUNDANCE_FILES

        for index,input_matrix in enumerate(input_matrices) : 

            input_matrix = os.path.join(options.features_dir,input_matrix)
            count_type,abundance_type = input_matrix.split('_')[1:3]
            check_file_exists(input_matrix)
            counts_all = {}
            header = []
            
            with open(input_matrix) as f:
                for line in f:
                    line = line.rstrip()
                    line_list = line.split('\t')
                    if len(header) == 0 :
                        header=line_list
                        for i in range(3,len(header),1) :
                            sample = header[i]
                            for annotation_id in annotation_id_list :
                                counts[annotation_id][sample] = 0
                            counts_all[sample] = 0 

                    else :
                        features = line_list[0]
                        annotation_id = features2annotation[features]
                        if annotation_id not in counts :
                            if annotation_id not in missing :
                                self.logger.warning("'%s' not present in %s" % (annotation_id,options.annotation_description))
                                missing.append(annotation_id)
                            continue
                        for i in range(3,len(header),1) :
                            sample = header[i]
                            counts[annotation_id][sample] = counts[annotation_id][sample] + float(line_list[i])
                            counts_all[sample] = counts_all[sample] + float(line_list[i])
            
            output_matrix = os.path.join(options.features_dir,output_matrices[index])
            self.logger.info('Print %s %s abundance matrix in "%s"' % (count_type, abundance_type, output_matrix))
            output_handle = open(output_matrix, "w")
            output_handle.write('\t'.join(['Features'] + header[3:len(header)])+'\n')
            for annotation in annotation_id_list :
                if sum([counts[annotation][s] for s in counts[annotation]]) == 0 and options.removed :
                    continue
                else :
                    output_handle.write('\t'.join([annotation] + [str(counts[annotation][s]) for s in counts[annotation]]) + '\n' )           
        
        self.logger.info('Printing matrices done')

    def parse_options(self, options):
        """Parse user options and call the correct pipeline(s)"""

        if(options.subparser_name == 'features'):
            self.features(options)
        elif(options.subparser_name == 'annoted_features'):
            self.annoted_features(options)
        elif(options.subparser_name == 'mama_wf'):
            options.features_dir = options.output_dir

            self.features(options)
            self.annoted_features(options)
        else:
            self.logger.error('Unknown MAMa command: ' + options.subparser_name + '\n')
            sys.exit(1)

        return 0
