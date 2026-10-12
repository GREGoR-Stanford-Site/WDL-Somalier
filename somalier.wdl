version 1.0

workflow Somalier {
  Array[Map[String,String]] omeList
  File DNA_sites
  File RNA_sites
  File reference
  File? pedigree

 scatter (ome in omeList) {
    Array[Array[File]] table = read_tsv(ome["toExtractList"])
    Array[File]? pastExtracted = if defined(ome["extractedList"]) then read_lines(ome["extractedList"]) else []
    
    # Determine which sites file to use based on ome name
    File sites_to_use = if (sub(ome["ome_name"], ".*RNA.*", "RNA") == "RNA") then RNA_sites else DNA_sites

    Float min_AB_threshold = if (sub(ome["ome_name"], ".*RNA.*", "RNA") == "RNA") then 0.2 else 0.3
    
    scatter (row in table) {
        call ExtractSample { 
            input: 
                sampleId=row[0], 
                sites=sites_to_use, 
                reference=reference,
                sampleBam=row[1], 
                sampleIndex=row[2]
        }
    }
    call RelateSamples {
        input:
            extractedFiles=ExtractSample.extractedFiles,
            oldExtractedFiles=pastExtracted,
            ome=ome["ome_name"],
            pedigree=pedigree,
            min_AB=min_AB_threshold
    }
    call CheckIdentical {
        input:
            pairsFile=RelateSamples.somalier_pairs,
            identityThreshold=0.95
    }
    if(defined(pedigree)) {
        call CheckRelationships {
            input:
                pairsFile=RelateSamples.somalier_pairs,
                relatednessThreshold=0.1
        }

        call CheckSex {
            input:
                sampleFile=RelateSamples.somalier_samples
        }
    }
  }

  # Cross-ome analysis: re-extract all samples using RNA sites for consistency
  Array[Array[Array[File]]] all_tables = table
  Array[Array[File]] flattened_bam_data = flatten(all_tables)
  
  scatter (sample_data in flattened_bam_data) {
    call ExtractSampleForCrossOme {
      input:
        sampleId=sample_data[0],
        sites=RNA_sites,  # Use RNA sites consistently for all cross-ome analysis
        reference=reference,
        sampleBam=sample_data[1],
        sampleIndex=sample_data[2]
    }
  }
  
  call CrossOmeAnalysis {
    input:
      extractedFiles=ExtractSampleForCrossOme.extractedFiles,
      pedigree=pedigree,
      min_AB=0.2  # Use RNA threshold for cross-ome analysis
  }
  
  call CheckCrossOmeSwaps {
    input:
      pairsFile=CrossOmeAnalysis.somalier_pairs,
      identityThreshold=0.95
  }

}

task ExtractSample {
  input {
    String sampleId
    File sites
    File reference
    File sampleBam
    File sampleIndex
  }
  command {
    SOMALIER_SAMPLE_NAME=${sampleId} somalier extract --sites ${sites} --fasta ${reference} ${sampleBam} 
  }
  output {
    # Write output to standard out
     File extractedFiles = "~{sampleId}.somalier"

  }
  runtime {
    docker: "quay.io/biocontainers/somalier:0.2.19--h0c29559_0"
  }
}

task ExtractSampleForCrossOme {
  input {
    String sampleId
    File sites
    File reference
    File sampleBam
    File sampleIndex
  }
  command {
    SOMALIER_SAMPLE_NAME=${sampleId}_cross_ome somalier extract --sites ${sites} --fasta ${reference} ${sampleBam} 
  }
  output {
    # Write output with cross_ome suffix to avoid name conflicts
     File extractedFiles = "~{sampleId}_cross_ome.somalier"

  }
  runtime {
    docker: "quay.io/biocontainers/somalier:0.2.19--h0c29559_0"
  }
}

task RelateSamples {
    Array[File] extractedFiles
    Array[File]? oldExtractedFiles
    File? pedigree
    String ome
    Float min_AB # minimum allele balance
    command {
        somalier relate -o ${ome} --min-ab ${min_AB} ${if defined(pedigree) then "-p ${pedigree}" else ""} ${sep=" " extractedFiles} ${if defined(oldExtractedFiles) then "${sep=' ' oldExtractedFiles}" else ""}
    }
    runtime {
        docker: "quay.io/biocontainers/somalier:0.2.19--h0c29559_0"
    }
    output {
        File somalier_pairs = "~{ome}.pairs.tsv"
        File somalier_samples = "~{ome}.samples.tsv"
        File somalier_html = "~{ome}.html"
    }
}

task CheckIdentical {
    File pairsFile
    Float identityThreshold
    command {
        python -c "
        import pandas as pd
        df = pd.read_csv('${pairsFile}',sep='\t')
        df = df[df['#sample_a'] == df['sample_b']]
        print('Sample Relatedness')
        for index, row in df.iterrows():
            if row['relatedness'] < ${identityThreshold}:
                print(row['#sample_a']+' '+str(row['relatedness']))
        "
    }
    runtime {
        docker: "quay.io/biocontainers/pandas:1.5.2"
    }
    output {
        File IdentityCheck = stdout()
    }
}

task CheckSex {
    File sampleFile
    command {
        python -c "
        import pandas as pd
        df = pd.read_csv('${sampleFile}',sep='\t')
        print('SampleId PredictedSex PedigreeSex')
        for index, row in df.iterrows():
            if (row['sex']==1 and row['original_pedigree_sex']=='male') or (row['sex']==2 and row['original_pedigree_sex']=='female'):
                continue
            else:
                print(row['sample_id'] + ' ' + str(row['sex']) + ' ' + row['original_pedigree_sex'])
        "
    }
    runtime {
        docker: "quay.io/biocontainers/pandas:1.5.2"
    }
    output {
        File SexCheck = stdout()
    }
}

task CheckRelationships {
    File pairsFile
    Float relatednessThreshold
    command {
        python -c "
        import pandas as pd 
        df = pd.read_csv('${pairsFile}',sep='\t')
        df = df[~(df['#sample_a'] == df['sample_b'])]
        print('Pair Relatedness Expected')
        for index, row in df.iterrows():
            if row['relatedness'] > row['expected_relatedness'] + ${relatednessThreshold} or row['relatedness'] < row['expected_relatedness'] - ${relatednessThreshold}:
                print(row['#sample_a']+'-'+row['sample_b']+' '+str(row['relatedness']) + ' '+ str(row['expected_relatedness']))
        "
    }
    runtime {
        docker: "quay.io/biocontainers/pandas:1.5.2"
    }
    output {
        File RelationshipCheck = stdout()
    }
}

task CrossOmeAnalysis {
    Array[File] extractedFiles
    File? pedigree
    Float min_AB
    command {
        somalier relate -o cross_ome_analysis --min-ab ${min_AB} ${sep=" " extractedFiles}
    }
    runtime {
        docker: "quay.io/biocontainers/somalier:0.2.19--h0c29559_0"
    }
    output {
        File somalier_pairs = "cross_ome_analysis.pairs.tsv"
        File somalier_samples = "cross_ome_analysis.samples.tsv"
        File somalier_html = "cross_ome_analysis.html"
    }
}

task CheckCrossOmeSwaps {
    File pairsFile
    Float identityThreshold
    command {
        python -c "
        import pandas as pd
        import re
        
        df = pd.read_csv('${pairsFile}', sep='\t')
        threshold = ${identityThreshold}
        
        # Extract base sample IDs (remove ome suffixes/prefixes)
        def extract_sample_id(sample_name):
            # Remove common ome prefixes/suffixes - adjust regex as needed for your naming convention
            base_id = re.sub(r'_(RNA|DNA|WGS|WES).*\$', '', sample_name)
            base_id = re.sub(r'^(RNA|DNA|WGS|WES)_', '', base_id)
            return base_id
        
        df['base_sample_a'] = df['#sample_a'].apply(extract_sample_id)
        df['base_sample_b'] = df['sample_b'].apply(extract_sample_id)
        
        print('Cross-Ome Sample Swap Analysis')
        print('==============================')
        print('Expected same individual, different omes (should have high relatedness):')
        
        # Check samples that should be the same individual across different omes
        same_individual = df[(df['base_sample_a'] == df['base_sample_b']) & 
                           (df['#sample_a'] != df['sample_b'])]
        
        for index, row in same_individual.iterrows():
            if row['relatedness'] < threshold:
                print('POTENTIAL SWAP: ' + row['#sample_a'] + ' vs ' + row['sample_b'] + ' - Relatedness: ' + str(round(row['relatedness'], 3)) + ' (Expected: >' + str(threshold) + ')')
        
        print()
        print('High relatedness between different individuals (potential swaps):')
        
        # Check for unexpectedly high relatedness between different individuals
        different_individuals = df[(df['base_sample_a'] != df['base_sample_b'])]
        
        for index, row in different_individuals.iterrows():
            if row['relatedness'] > threshold:
                print('POTENTIAL SWAP: ' + row['#sample_a'] + ' vs ' + row['sample_b'] + ' - Relatedness: ' + str(round(row['relatedness'], 3)) + ' (Expected: <' + str(threshold) + ')')
        "
    }
    runtime {
        docker: "quay.io/biocontainers/pandas:1.5.2"
    }
    output {
        File CrossOmeSwapCheck = stdout()
    }
}
