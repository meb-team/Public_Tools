#########################################################################################
#                                                                                       #
# defaultValues.py - store default values used in many places in CheckM                 #
#                                                                                       #
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

class DefaultValues():
    """Default values for filenames and common constants."""

    FEATURES_ABUNDANCE_FILES = ['features_reads_raw_abundance.tsv','features_reads_normalised_abundance.tsv','features_reads_relative_abundance.tsv','features_base_raw_abundance.tsv','features_base_normalised_abundance.tsv','features_base_relative_abundance.tsv' ]

    ANNOTATE_ABUNDANCE_FILES = ['annotate_reads_raw_abundance.tsv', 'annotate_reads_normalised_abundance.tsv', 'annotate_reads_relative_abundance.tsv', 'annotate_base_raw_abundance.tsv', 'annotate_base_normalised_abundance.tsv' , 'annotate_base_relative_abundance.tsv']

