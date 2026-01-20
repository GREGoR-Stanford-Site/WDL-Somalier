version 1.0

workflow Somalier {
  input {
    Array[Map[String,String]] omeList
    File sites
    File reference
    File? pedigree
  }

  scatter (ome in omeList) {
    Array[Array[File]] table = read_tsv(ome["toExtractList"])
    Array[File] pastExtracted = read_lines(ome["extractedList"])
    scatter (row in table) {
        call ExtractSample { 
            input: 
                sampleId=row[0], 
                sites=sites, 
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
            pedigree=pedigree
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

  # Check identity across all omes
  call CheckIdenticalAcrossOmes {
    input:
      pairsFiles=RelateSamples.somalier_pairs,
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

task RelateSamples {
    input {
        Array[File] extractedFiles
        Array[File]? oldExtractedFiles
        File? pedigree
        String ome
    }
    command {
        somalier relate -o ${ome} ${if defined(pedigree) then "-p ${pedigree}" else ""} ${sep=" " extractedFiles} ${if defined(oldExtractedFiles) then "sep=' ' oldExtractedFiles" else ""}
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
    input {
        File pairsFile
        Float identityThreshold
    }
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
    input {
        File sampleFile
    }
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
    input {
        File pairsFile
        Float relatednessThreshold
    }
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

task CheckIdenticalAcrossOmes {
    input {
        Array[File] pairsFiles
        Float identityThreshold
    }
    command <<<
        python -c "
import pandas as pd
import sys
# Combine all pairs files from different omes
dfs = []
for f in sys.argv[1:]:
    dfs.append(pd.read_csv(f, sep='\t'))
if not dfs:
    print('Sample Relatedness (Across All Omes)')
    sys.exit(0)
df = pd.concat(dfs, ignore_index=True)
# Filter for self-comparisons (same sample)
df = df[df['#sample_a'] == df['sample_b']]
print('Sample Relatedness (Across All Omes)')
for index, row in df.iterrows():
    if row['relatedness'] < ~{identityThreshold}:
        print(row['#sample_a']+' '+str(row['relatedness']))
" ~{sep=" " pairsFiles}
    >>>
    runtime {
        docker: "quay.io/biocontainers/pandas:1.5.2"
    }
    output {
        File IdentityCheckAcrossOmes = stdout()
    }
}
