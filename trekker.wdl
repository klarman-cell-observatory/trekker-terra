version 1.0

workflow trekker_sc {
    input {
        File input_samplesheet
        String output_directory
        String docker_registry
        Boolean run_trekker_merge
        Boolean existing_trekker_outs

        Int num_cpu = 16
        String memory = "256G"
        Int disk_space = 300
    }

    Array[String] all_rows = read_lines(input_samplesheet)

    if (!existing_trekker_outs) {
        scatter (row_index in range(length(all_rows) - 1)) {
            Int row_index_plus_one = row_index + 1
            String sample_line = all_rows[row_index_plus_one]

            call process_sample {
                input:
                    line = sample_line,
                    input_samplesheet = input_samplesheet,
                    output_directory = output_directory,
                    docker_registry = docker_registry,
                    num_cpu = num_cpu,
                    memory = memory,
                    disk_space = disk_space
            }
        }

        if (run_trekker_merge) {
            call merge_sample_groups as merge_after_run {
                input:
                    input_samplesheet = input_samplesheet,
                    output_directory = output_directory,
                    docker_registry = docker_registry,
                    run_outputs = flatten(select_first([process_sample.trekker_output_files, []])),
                    existing_trekker_outs = false,
                    num_cpu = num_cpu,
                    memory = memory,
                    disk_space = disk_space
            }
        }
    }

    if (existing_trekker_outs && run_trekker_merge) {
        call merge_sample_groups as merge_existing {
            input:
                input_samplesheet = input_samplesheet,
                output_directory = output_directory,
                docker_registry = docker_registry,
                run_outputs = [],
                existing_trekker_outs = true,
                num_cpu = num_cpu,
                memory = memory,
                disk_space = disk_space
        }
    }

    output {
        Array[File] trekker_outputs = if (existing_trekker_outs) then [] else flatten(select_first([process_sample.trekker_output_files, []]))
    }
}


task process_sample {
    input {
        String line
        File input_samplesheet
        String output_directory
        String docker_registry
        Int num_cpu
        String memory
        Int disk_space
    }

    command <<<
        set -euo pipefail

        MNT_PATH="/mnt/disks/cromwell_root"
        RAW_DIR="${MNT_PATH}/raw"
        OUT_DIR="${MNT_PATH}/out"
        mkdir -p "${RAW_DIR}" "${OUT_DIR}"
        chmod -R 777 "${MNT_PATH}"

        python3 - \
            "~{input_samplesheet}" \
            "~{line}" \
            "${RAW_DIR}/sample_info.tsv" <<'PY'
import csv, sys
sheet, line, output = sys.argv[1:]
with open(sheet, newline='') as f:
    reader = csv.DictReader(f)
    header = reader.fieldnames or []
    if not header:
        raise RuntimeError("Samplesheet is empty or has no header")
parsed = next(csv.reader([line]))
if len(parsed) != len(header):
    raise RuntimeError(f"Samplesheet row has {len(parsed)} columns; expected {len(header)}")
data = dict(zip(header, parsed))
required = ["sample", "sample_group", "barcode_file", "fastq_1", "fastq_2", "sc_outdir"]
missing = [k for k in required if k not in data]
if missing:
    raise RuntimeError("Samplesheet missing required columns: " + ", ".join(missing))
with open(output, "w") as f:
    for key in required:
        f.write(f"{key}\t{data[key]}\n")
PY

        get_value() { awk -F '\t' -v key="$1" '$1 == key {sub(/^[^\t]*\t/, ""); print; exit}' "${RAW_DIR}/sample_info.tsv"; }
        SAMPLE_NAME="$(get_value sample)"
        SAMPLE_GROUP="$(get_value sample_group)"
        BARCODE_PATH="$(get_value barcode_file)"
        FASTQR1_PATH="$(get_value fastq_1)"
        FASTQR2_PATH="$(get_value fastq_2)"
        SC_OUTDIR="$(get_value sc_outdir)"

        if [[ -z "${SAMPLE_NAME}" || -z "${SAMPLE_GROUP}" ]]; then
            echo "ERROR: sample and sample_group must both be populated in the samplesheet." >&2
            exit 1
        fi

        echo "Processing sample: ${SAMPLE_NAME} (group: ${SAMPLE_GROUP})"
        echo "R1: ${FASTQR1_PATH}"
        echo "R2: ${FASTQR2_PATH}"
        echo "Barcode: ${BARCODE_PATH}"
        echo "sc_outdir: ${SC_OUTDIR}"

        # Localize the sample's inputs.
        gcloud storage cp "${FASTQR1_PATH}" "${RAW_DIR}/"
        gcloud storage cp "${FASTQR2_PATH}" "${RAW_DIR}/"
        gcloud storage cp "${BARCODE_PATH}" "${RAW_DIR}/"
        LOCAL_R1="${RAW_DIR}/$(basename "${FASTQR1_PATH}")"
        LOCAL_R2="${RAW_DIR}/$(basename "${FASTQR2_PATH}")"
        LOCAL_BARCODE="${RAW_DIR}/$(basename "${BARCODE_PATH}")"

        if [[ "${SC_OUTDIR}" == *.h5 ]]; then
            gcloud storage cp "${SC_OUTDIR}" "${RAW_DIR}/"
            LOCAL_CELLBENDER_H5="${RAW_DIR}/$(basename "${SC_OUTDIR}")"
            LOCAL_SC_OUTDIR="${RAW_DIR}/cellbender_matrix"
            python3 /opt/trekker-v1.4.11/terra/cellbender_to_10x.py \
                "${LOCAL_CELLBENDER_H5}" "${LOCAL_SC_OUTDIR}"
        else
            LOCAL_SC_PARENT="${RAW_DIR}/sc_input"
            mkdir -p "${LOCAL_SC_PARENT}"
            gcloud storage cp -r "${SC_OUTDIR}" "${LOCAL_SC_PARENT}/"
            LOCAL_SC_OUTDIR="${LOCAL_SC_PARENT}/$(basename "${SC_OUTDIR%/}")"
        fi

        # Write a one-sample samplesheet with localized paths.
        LOCAL_TREKKER_SAMPLESHEET="${RAW_DIR}/trekker_samplesheet.csv"
        python3 - \
            "~{input_samplesheet}" \
            "${LOCAL_TREKKER_SAMPLESHEET}" \
            "${SAMPLE_NAME}" \
            "${LOCAL_BARCODE}" \
            "${LOCAL_R1}" \
            "${LOCAL_R2}" \
            "${LOCAL_SC_OUTDIR}" <<'PY'
import csv, sys
source, destination, sample, barcode, r1, r2, sc_outdir = sys.argv[1:]
with open(source, newline='') as f:
    reader = csv.DictReader(f)
    if not reader.fieldnames:
        raise RuntimeError("Samplesheet has no header")
    header = reader.fieldnames
    required = ["sample", "barcode_file", "fastq_1", "fastq_2", "sc_outdir"]
    missing = [c for c in required if c not in header]
    if missing:
        raise RuntimeError("Samplesheet missing required columns: " + ", ".join(missing))
    rows = list(reader)
matched = [r for r in rows if r.get("sample") == sample]
if len(matched) != 1:
    raise RuntimeError(f"Expected exactly one row for sample {sample}; found {len(matched)}")
row = matched[0]
row["barcode_file"] = barcode
row["fastq_1"] = r1
row["fastq_2"] = r2
row["sc_outdir"] = sc_outdir
with open(destination, "w", newline='') as f:
    writer = csv.DictWriter(f, fieldnames=header)
    writer.writeheader()
    writer.writerow(row)
print(f"Wrote localized samplesheet: {destination}")
PY

        cat "${LOCAL_TREKKER_SAMPLESHEET}"
        bash /opt/trekker-v1.4.11/nuclei_locater_toplevel.sh "${LOCAL_TREKKER_SAMPLESHEET}"

        ANALYSIS_LOG="${RAW_DIR}/log/${SAMPLE_NAME}/analysis.log"
        if [[ -f "${ANALYSIS_LOG}" ]]; then
            cp "${ANALYSIS_LOG}" "${OUT_DIR}/${SAMPLE_NAME}_analysis.log"
        fi

        OUTPUT_FILE_COUNT="$(find "${OUT_DIR}" -type f | wc -l)"
        if [[ "${OUTPUT_FILE_COUNT}" -eq 0 ]]; then
            echo "ERROR: Trekker produced no files for ${SAMPLE_NAME}." >&2
            exit 1
        fi

        DESTINATION="~{output_directory}"
        DESTINATION="${DESTINATION%/}/${SAMPLE_NAME}/outputs"
        while IFS= read -r FILE; do
            RELATIVE_PATH="${FILE#${OUT_DIR}/}"
            gcloud storage cp "${FILE}" "${DESTINATION}/${RELATIVE_PATH}"
        done < <(find "${OUT_DIR}" -type f -print)

        echo "Uploaded ${OUTPUT_FILE_COUNT} output files to ${DESTINATION}"
    >>>

    output {
        Array[File] trekker_output_files = glob("/mnt/disks/cromwell_root/out/**/*")
    }

    runtime {
        docker: docker_registry
        cpu: num_cpu
        memory: memory
        bootDiskSizeGB: 25
        disks: "local-disk ${disk_space} HDD"
    }
}


task merge_sample_groups {
    input {
        File input_samplesheet
        String output_directory
        String docker_registry
        Array[File] run_outputs
        Boolean existing_trekker_outs
        Int num_cpu
        String memory
        Int disk_space
    }

    command <<<
        set -euo pipefail

        ROOT="/mnt/disks/cromwell_root"
        WORK_DIR="${ROOT}/merge_work"
        mkdir -p "${WORK_DIR}" "${WORK_DIR}/merged_results"
        chmod -R 777 "${ROOT}"

        # Read the samplesheet and create a manifest of sample/group pairs.
        # The merger samplesheet header is sample,out_dir.
        python3 - "~{input_samplesheet}" "${WORK_DIR}/sample_manifest.tsv" <<'PY'
import csv, sys
source, destination = sys.argv[1:]
with open(source, newline='') as f:
    reader = csv.DictReader(f)
    if not reader.fieldnames:
        raise RuntimeError("Samplesheet is empty or has no header")
    required = ["sample", "sample_group"]
    missing = [c for c in required if c not in reader.fieldnames]
    if missing:
        raise RuntimeError("Samplesheet missing required columns: " + ", ".join(missing))
    rows = list(reader)
if not rows:
    raise RuntimeError("Samplesheet contains no sample rows")
seen = set()
with open(destination, "w") as out:
    for row in rows:
        sample = (row.get("sample") or "").strip()
        group = (row.get("sample_group") or "").strip()
        if not sample or not group:
            raise RuntimeError("Every row must have non-empty sample and sample_group values")
        if (sample, group) in seen:
            raise RuntimeError(f"Duplicate sample/group row: {sample}/{group}")
        seen.add((sample, group))
        out.write(sample + "\t" + group + "\n")
PY

        OUTPUT_BASE="~{output_directory}"
        OUTPUT_BASE="${OUTPUT_BASE%/}"
        MANIFEST="${WORK_DIR}/sample_manifest.tsv"
        MERGE_SCRIPT="/opt/trekker_merger-v1.4.11/trekker_merger.sh"

        if [[ ! -x "${MERGE_SCRIPT}" && ! -f "${MERGE_SCRIPT}" ]]; then
            echo "ERROR: expected merger script at ${MERGE_SCRIPT}. Ensure it is included in the Docker image." >&2
            exit 1
        fi

        cut -f2 "${MANIFEST}" | awk '!seen[$0]++' > "${WORK_DIR}/groups.txt"
        : > "${WORK_DIR}/merge_manifest.txt"

        while IFS= read -r GROUP; do
            [[ -n "${GROUP}" ]] || continue

            SAFE_GROUP="$(printf '%s' "${GROUP}" | sed 's/[^a-zA-Z0-9_]/_/g')"
            GROUP_DIR="${WORK_DIR}/groups/${SAFE_GROUP}"
            mkdir -p "${GROUP_DIR}/samples" "${GROUP_DIR}/merged_results"

            MERGE_SAMPLESHEET="${GROUP_DIR}/merge_samplesheet.csv"
            printf 'sample,out_dir\n' > "${MERGE_SAMPLESHEET}"

            while IFS=$'\t' read -r SAMPLE SAMPLE_GROUP; do
                [[ "${SAMPLE_GROUP}" == "${GROUP}" ]] || continue

                LOCAL_SAMPLE_PARENT="${GROUP_DIR}/samples/${SAMPLE}"
                mkdir -p "${LOCAL_SAMPLE_PARENT}"

                if [[ "~{existing_trekker_outs}" == "true" ]]; then
                    # Search for one known file, rather than listing the entire
                    # output directory. Derive the output directory from that file.
                    SAMPLE_ROOT="${OUTPUT_BASE}/${SAMPLE}"
                    MATCH_LIST="${GROUP_DIR}/${SAMPLE}_metrics_matches.txt"
                    METRICS_NAME="${SAMPLE}_summary_metrics.csv"
                    RDS_NAME="${SAMPLE}_seurat_spatial.rds"

                    gcloud storage ls \
                        "${SAMPLE_ROOT}/*/trekker_${SAMPLE}/output/${METRICS_NAME}" \
                        > "${MATCH_LIST}" || true

                    MATCH_COUNT="$(grep -c . "${MATCH_LIST}" || true)"
                    if [[ "${MATCH_COUNT}" -ne 1 ]]; then
                        echo "ERROR: expected exactly one existing summary metrics file for ${SAMPLE}; found ${MATCH_COUNT}." >&2
                        echo "Search pattern: ${SAMPLE_ROOT}/*/trekker_${SAMPLE}/output/${METRICS_NAME}" >&2
                        cat "${MATCH_LIST}" >&2 || true
                        exit 1
                    fi

                    CLOUD_METRICS_FILE="$(head -n 1 "${MATCH_LIST}" | sed 's/[[:space:]]*$//')"
                    CLOUD_OUTPUT_DIR="$(dirname "${CLOUD_METRICS_FILE}")"
                    CLOUD_RDS_FILE="${CLOUD_OUTPUT_DIR}/intermediates/${RDS_NAME}"

                    # Recreate only the directory structure required by the merger.
                    LOCAL_SAMPLE_OUT="${LOCAL_SAMPLE_PARENT}/output"
                    mkdir -p "${LOCAL_SAMPLE_OUT}/intermediates"

                    echo "Localizing required merger files for ${SAMPLE} from ${CLOUD_OUTPUT_DIR}"
                    gcloud storage cp \
                        "${CLOUD_METRICS_FILE}" \
                        "${LOCAL_SAMPLE_OUT}/${METRICS_NAME}"

                    gcloud storage cp \
                        "${CLOUD_RDS_FILE}" \
                        "${LOCAL_SAMPLE_OUT}/intermediates/${RDS_NAME}"

                else
                    # In run mode, use the complete output tree produced by this workflow.
                    CLOUD_SAMPLE_DIR="${OUTPUT_BASE}/${SAMPLE}/outputs/"
                    LOCAL_SAMPLE_OUT="${LOCAL_SAMPLE_PARENT}/outputs"

                    echo "Localizing ${CLOUD_SAMPLE_DIR} to ${LOCAL_SAMPLE_OUT}"
                    gcloud storage cp -r "${CLOUD_SAMPLE_DIR}" "${LOCAL_SAMPLE_PARENT}/"

                    if [[ ! -d "${LOCAL_SAMPLE_OUT}" ]]; then
                        echo "ERROR: expected localized sample output directory ${LOCAL_SAMPLE_OUT}" >&2
                        find "${LOCAL_SAMPLE_PARENT}" -maxdepth 5 -print
                        exit 1
                    fi
                fi

                printf '%s,%s\n' "${SAMPLE}" "${LOCAL_SAMPLE_OUT}" >> "${MERGE_SAMPLESHEET}"

            done < "${MANIFEST}"

            MERGED_NAME="${SAFE_GROUP}_merged"
            MERGED_LOCAL="${GROUP_DIR}/merged_results/${MERGED_NAME}"
            mkdir -p "$(dirname "${MERGED_LOCAL}")"

            echo "Merging group '${GROUP}' into ${MERGED_NAME}"
            bash "${MERGE_SCRIPT}" \
                "${MERGE_SAMPLESHEET}" \
                "${MERGED_LOCAL}" \
                "${MERGED_NAME}" \
                "conda"

            gcloud storage cp -r "${MERGED_LOCAL}" "${OUTPUT_BASE}/"
            echo "${GROUP}	${OUTPUT_BASE}/${MERGED_NAME}" >> "${WORK_DIR}/merge_manifest.txt"

        done < "${WORK_DIR}/groups.txt"

        cp "${WORK_DIR}/merge_manifest.txt" "${WORK_DIR}/merge_manifest.txt.copy"
        echo "Merge complete. Group output locations:"
        cat "${WORK_DIR}/merge_manifest.txt"
    >>>

    output {
        File merge_manifest = "/mnt/disks/cromwell_root/merge_work/merge_manifest.txt.copy"
    }

    runtime {
        docker: docker_registry
        cpu: num_cpu
        memory: memory
        bootDiskSizeGB: 25
        disks: "local-disk ${disk_space} HDD"
    }
}