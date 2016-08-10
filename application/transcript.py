#FASTA Datapull to memory
#Includes transcript and 5000bp upstream and downstream
def ensembl_sequence(ensembl_transcript_id):
    import requests
    query_in_string = \
        '''http://ensembl.org/biomart/martservice?query='''\
        '''<?xml version="1.0" encoding="UTF-8"?>'''\
        '''<!DOCTYPE Query>'''\
        '''<Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >'''\
        '''<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'''\
        '''<Filter name = "ensembl_transcript_id" value = "%s"/>'''\
        '''<Filter name = "downstream_flank" value = "5000"/>'''\
        '''<Filter name = "upstream_flank" value = "5000"/>'''\
        '''<Attribute name = "ensembl_transcript_id" />'''\
        '''<Attribute name = "transcript_exon_intron" />'''\
        '''</Dataset>'''\
        '''</Query>'''\

    text = 'Q'


    #HTTP Query loop for the ensembl_flank_sequence filter not found error. This has to be reworked to prevent HTTP blocking
    while text[0] == 'Q':
        url = query_in_string % (ensembl_transcript_id)
        req = requests.get(url, stream=True)
        text = req.text
    return text.lower()


#Ensembl table pull: Transcript Summary, SNP table, Exon/UTR table
def ensembl_summary(ensembl_transcript_id):
    import requests
    import StringIO
    import csv

    query_in_string = \
    '''http://ensembl.org/biomart/martservice?query='''\
    '''<?xml version="1.0" encoding="UTF-8"?>'''\
    '''<!DOCTYPE Query>'''\
    '''<Query  virtualSchemaName = "default" formatter = "CSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >'''\
    '''<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'''\
    '''<Filter name = "ensembl_transcript_id" value = "%s"/>'''\
    '''<Attribute name = "ensembl_transcript_id" />'''\
    '''<Attribute name = "external_transcript_name" />'''\
    '''<Attribute name = "chromosome_name" />'''\
    '''<Attribute name = "transcript_start" />'''\
    '''<Attribute name = "transcript_end" />'''\
    '''<Attribute name = "strand" />'''\
    '''<Attribute name = "transcript_biotype" />'''\
    '''<Attribute name = "transcript_source" />'''\
    '''<Attribute name = "transcript_status" />'''\
    '''<Attribute name = "transcript_version" />'''\
    '''</Dataset>'''\
    '''</Query>'''\

    url = query_in_string % (ensembl_transcript_id)
    req = requests.get(url, stream=True)
    #print req.text
    input = StringIO.StringIO(req.text)
    reader = csv.DictReader(input,delimiter = ',')
    list = []
    for row in reader:
        list.append(row)
    return list

def ensembl_snp(ensembl_transcript_id):
    import requests
    import StringIO
    import csv

    query_in_string = \
    '''http://ensembl.org/biomart/martservice?query='''\
    '''<?xml version="1.0" encoding="UTF-8"?>'''\
    '''<!DOCTYPE Query>'''\
    '''<Query  virtualSchemaName = "default" formatter = "CSV" header = "1" uniqueRows = "1" count = "" datasetConfigVersion = "0.6" >'''\
    '''<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'''\
    '''<Filter name = "ensembl_transcript_id" value = "%s"/>'''\
    '''<Attribute name = "variation_name" />'''\
    '''<Attribute name = "chromosome_name" />'''\
    '''<Attribute name = "chromosome_start" />'''\
    '''<Attribute name = "chromosome_end" />'''\
    '''<Attribute name = "strand" />'''\
    '''<Attribute name = "synonymous_status" />'''\
    '''<Attribute name = "allele" />'''\
    '''<Attribute name = "minor_allele_freq" />'''\
    '''</Dataset>'''\
    '''</Query>'''\

    url = query_in_string % (ensembl_transcript_id)
    req = requests.get(url, stream=True)
    input = StringIO.StringIO(req.text)
    reader = csv.DictReader(input,delimiter = ',')
    list = []
    for row in reader:
        list.append(row)
    return list

def ensembl_exon_utr(ensembl_transcript_id):
    import requests
    import StringIO
    import csv


    query_in_string = \
    '''http://ensembl.org/biomart/martservice?query='''\
    '''<?xml version="1.0" encoding="UTF-8"?>'''\
    '''<!DOCTYPE Query>'''\
    '''<Query  virtualSchemaName = "default" formatter = "CSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >'''\
    '''<Dataset name = "hsapiens_gene_ensembl" interface = "default" >'''\
    '''<Filter name = "ensembl_transcript_id" value = "%s"/>'''\
    '''<Attribute name = "exon_chrom_start" />'''\
    '''<Attribute name = "exon_chrom_end" />'''\
    '''<Attribute name = "3_utr_start" />'''\
    '''<Attribute name = "3_utr_end" />'''\
    '''<Attribute name = "5_utr_start" />'''\
    '''<Attribute name = "5_utr_end" />'''\
    '''</Dataset>'''\
    '''</Query>'''\

    url = query_in_string % (ensembl_transcript_id)
    req = requests.get(url, stream=True)
    input = StringIO.StringIO(req.text)
    reader = csv.DictReader(input,delimiter = ',')
    list = []
    for row in reader:
        list.append(row)
    return list

#Transcript Highlight Methods, forward and reverse
def highlight(sequence,summary,exon_utr,snp,bisulfite):
    import re
    from bs4 import BeautifulSoup

    #Sequence prep
    sequence = sequence[17:-1].replace('\n','').replace('\t','')
    sequence = list(sequence)

    #sequence snp and UTR highlights

    start_utr_tag = '''<span style="color:brown;">'''
    end_utr_tag = "</span>"

    start_snp_link = '''<a href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=%s">'''
    end_snp_link = "</a>"

    start_snp_tag_utr = '''<span style="background-color:#7ac5cd;color:#000000;">'''
    start_snp_tag_downstream = '''<span style="background-color:#a2b5cd;color:#000000;">'''
    start_snp_tag_frameshift= '''<span style="background-color:#9400D3;color:#ffffff;">'''
    start_snp_tag_intronic = '''<span style="background-color:#02599c;color:#ffffff;">'''
    start_snp_tag_missense = '''<span style="background-color:#ffd700;color:#000000;">'''
    start_snp_tag_splice_acceptor= '''<span style="background-color:#FF581A;color:#ffffff;">'''
    start_snp_tag_splice_region = '''<span style="background-color:#ff7f50;color:#ffffff;">'''
    start_snp_tag_stop_gained = '''<span style="background-color:#ff0000;color:#ffffff;">'''
    start_snp_tag_synonymous = '''<span style="background-color:#76ee00;color:#000000;">'''
    end_snp_tag = "</span>"

    start_cg_tag = '''<span style="color:red">'''
    end_cg_tag = "</span>"

    #CpG Site highlighting in full and bisulfite sequence
    for term in range(0,len(sequence)):
        if sequence[term] == 'c':
            if sequence[term+1] =='g':
                if bisulfite==False:
                    sequence[term] = start_cg_tag + sequence[term] + end_cg_tag
                    sequence[term+1] = start_cg_tag + sequence[term+1] + end_cg_tag
                else:
                    sequence[term] = start_cg_tag + 'y' +end_cg_tag
                    sequence[term+1] = start_cg_tag + sequence[term+1] + end_cg_tag

    #Sequence Exon capitalize; method changes depending on strand direction
    if int(summary[0]['Strand']) == -1:

        zero_index_offset = 0
        flank_offset = 5000
        end_term_offset = 1
        transcript_end = int(summary[0]['Transcript End (bp)'])
        total_offset = transcript_end + flank_offset + zero_index_offset

        for row in exon_utr:

            #Capitalize all Exons
            start_pos = total_offset - int(row['Exon Chr End (bp)'])
            end_pos = total_offset - int(row['Exon Chr Start (bp)'])
            for term in range(start_pos,end_pos+end_term_offset):
                sequence[term] = sequence[term].upper()


            #Color UTR Brown
            if row["5' UTR Start"].isdigit() == True:
                start_pos = total_offset - int(row["5' UTR End"])
                end_pos = total_offset - int(row["5' UTR Start"])
                for term in range(start_pos,end_pos):
                    sequence[term] = start_utr_tag + sequence[term] + end_utr_tag


            if row["3' UTR Start"].isdigit() == True:
                start_pos = total_offset - int(row["3' UTR End"])
                end_pos = total_offset - int(row["3' UTR Start"])
                for term in range(start_pos,end_pos):
                    sequence[term] = start_utr_tag + sequence[term] + end_utr_tag


        for row in snp:
            start_pos = total_offset - int(row['Chromosome position end (bp)'])
            end_pos = total_offset - int(row['Chromosome position start (bp)'])
            for term in range(start_pos,end_pos+end_term_offset):
                start_snp_link_updated = start_snp_link %(row['Variant Name'])
                sequence[term] = start_snp_link_updated + sequence[term]+ end_snp_link

                if row['Variant Consequence']=='downstream_gene_variant':
                    sequence[term] = start_snp_tag_downstream + sequence[term] + end_snp_tag

                elif row['Variant Consequence']=='upstream_gene_variant':
                    sequence[term] = start_snp_tag_downstream + sequence[term] + end_snp_tag

                elif row['Variant Consequence']=='3_prime_UTR_variant':
                    sequence[term] = start_snp_tag_utr + sequence[term] + end_snp_tag

                elif row['Variant Consequence']=='5_prime_UTR_variant':
                    sequence[term] = start_snp_tag_utr + sequence[term] + end_snp_tag

                elif row['Variant Consequence']=='synonymous_variant':
                    sequence[term] = start_snp_tag_synonymous + sequence[term] + end_snp_tag

                elif row['Variant Consequence']=='missense_variant':
                    sequence[term] = start_snp_tag_missense + sequence[term] + end_snp_tag

                elif row['Variant Consequence']=='splice_region_variant':
                    sequence[term] = start_snp_tag_splice_region + sequence[term] + end_snp_tag

                elif row['Variant Consequence']=='intron_variant':
                    sequence[term] = start_snp_tag_intronic + sequence[term] + end_snp_tag

                elif row['Variant Consequence']=='stop_gained':
                    sequence[term] = start_snp_tag_stop_gained + sequence[term] + end_snp_tag

                elif row['Variant Consequence']=='frameshift_variant':
                    sequence[term] = start_snp_tag_frameshift + sequence[term] + end_snp_tag

                elif row['Variant Consequence']=='splice_acceptor_variant':
                    sequence[term] = start_snp_tag_splice_acceptor + sequence[term] + end_snp_tag

                elif row['Variant Consequence']=='splice_donor_variant':
                    sequence[term] = start_snp_tag_splice_acceptor + sequence[term] + end_snp_tag




    #If forward strand
    else:
        zero_index_offset = 0
        flank_offset = 5000
        end_term_offset = 1
        transcript_start = int(summary[0]['Transcript Start (bp)'])
        total_offset = transcript_start - flank_offset + zero_index_offset

        for row in exon_utr:
            start_pos = int(row['Exon Chr Start (bp)']) - total_offset
            end_pos = int(row['Exon Chr End (bp)']) - total_offset
            for term in range(start_pos,end_pos+end_term_offset):
                sequence[term] = sequence[term].upper()


            #Color UTR Brown
            if row["5' UTR Start"].isdigit() == True:

                start_pos = int(row["5' UTR Start"]) - total_offset
                end_pos = int(row["5' UTR End"]) - total_offset
                for term in range(start_pos,end_pos):
                    sequence[term] = start_utr_tag + sequence[term] + end_utr_tag


            if row["3' UTR Start"].isdigit() == True:
                start_pos = int(row["3' UTR Start"]) - total_offset
                end_pos = int(row["3' UTR End"]) - total_offset
                for term in range(start_pos,end_pos):
                    sequence[term] = start_utr_tag + sequence[term] + end_utr_tag

        for row in snp:
            start_pos = int(row['Chromosome position end (bp)']) - total_offset
            end_pos = int(row['Chromosome position start (bp)']) - total_offset

            for term in range(start_pos,end_pos+end_term_offset):
                start_snp_link_updated = start_snp_link %(row['Variant Name'])
                sequence[term] = start_snp_link_updated + sequence[start_pos] + end_snp_link


                if row['Variant Consequence']=='downstream_gene_variant':
                    sequence[term] = start_snp_tag_downstream + sequence[term] + end_snp_tag

                elif row['Variant Consequence']=='upstream_gene_variant':
                    sequence[term] = start_snp_tag_downstream + sequence[term] + end_snp_tag

                elif row['Variant Consequence']=='3_prime_UTR_variant':
                    sequence[term] = start_snp_tag_utr + sequence[term] + end_snp_tag


                elif row['Variant Consequence']=='5_prime_UTR_variant':
                    sequence[term] = start_snp_tag_utr + sequence[term] + end_snp_tag


                elif row['Variant Consequence']=='synonymous_variant':
                    sequence[term] = start_snp_tag_synonymous + sequence[term] + end_snp_tag

                elif row['Variant Consequence']=='missense_variant':
                    sequence[term] = start_snp_tag_missense + sequence[term] + end_snp_tag

                elif row['Variant Consequence']=='splice_region_variant':
                    sequence[term] = start_snp_tag_splice_region + sequence[term] + end_snp_tag


                elif row['Variant Consequence']=='intron_variant':
                    sequence[term] = start_snp_tag_intronic + sequence[term] + end_snp_tag


                elif row['Variant Consequence']=='stop_gained':
                    sequence[term] = start_snp_tag_stop_gained + sequence[term] + end_snp_tag


                elif row['Variant Consequence']=='frameshift_variant':
                    sequence[term] = start_snp_tag_frameshift + sequence[term] + end_snp_tag


                elif row['Variant Consequence']=='splice_acceptor_variant':
                    sequence[term] = start_snp_tag_splice_acceptor + sequence[term] + end_snp_tag


                elif row['Variant Consequence']=='splice_donor_variant':
                    sequence[term] = start_snp_tag_splice_acceptor + sequence[term] + end_snp_tag



    soup =  BeautifulSoup(''.join(sequence),"lxml")
    full_sequence = soup.prettify()
    regex = re.compile(r"\s*(<[^<>]+>)\s*")
    full_sequence = regex.sub("\g<1>", full_sequence)
    return full_sequence



if __name__ == "__main__":
    ensembl_transcript_id = 'ENST00000217270'

    sequence = ensembl_sequence(ensembl_transcript_id)

    summary = ensembl_summary(ensembl_transcript_id)
    exon_utr = ensembl_exon_utr(ensembl_transcript_id)
    snp = ensembl_snp(ensembl_transcript_id)

    transcript = highlight(sequence,summary,exon_utr,snp,False)

    bisulfite = highlight(sequence,summary,exon_utr,snp,True)


    f = open('test.html',mode='w')
    f.write(transcript)
    f.close()


    g = open('test2.html',mode='w')
    g.write(bisulfite)
    g.close()
