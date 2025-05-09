#!/bin/bash
#
# 批量处理脚本，用于处理样本表中的样本并准备mageck输入文件
# 用途：批量处理sgRNA测序数据，并整理结果以符合mageck的输入要求
#
# 用法: $0 <samplesheet.csv> <crop_length> [选项]
#   samplesheet.csv - 样本表文件
#   crop_length    - 裁剪长度
#
# 选项:
#   --index=NAME   - 指定索引名称
#   --merge        - 合并gene和promoter数据
#   --run-mageck   - 自动运行mageck分析
#   --mode=MODE    - mageck分析模式 (test或mle，默认为test)
#

set -e  # 遇到错误立即退出
set -o pipefail  # 管道中任何命令失败都会导致整个管道失败

# 显示使用方法
function show_usage {
	echo "用法: $0 <samplesheet.csv> <crop_length> [选项]"
	echo "选项:"
	echo "  --index=NAME   - 指定索引名称"
	echo "  --merge        - 合并gene和promoter数据"
	echo "  --run-mageck   - 自动运行mageck分析"
	echo "  --mode=MODE    - mageck分析模式 (test或mle，默认为test)"
	echo "示例: $0 /home/gyk/project/SiCRIS/data/samplesheet.csv 100 --merge --run-mageck --mode=mle"
	exit 1
}

# 参数检查
if [ $# -lt 2 ]; then
	echo "错误: 参数不足"
	show_usage
fi

# 检查样本表文件是否存在
if [ ! -f "$1" ]; then
	echo "错误: 样本表文件 $1 不存在"
	exit 1
fi

SAMPLESHEET="$1"
CROP_LENGTH="$2"
INDEX_NAME=""
MERGE_INDICES=false
RUN_MAGECK=false
MAGECK_MODE="test"

# 解析可选参数
shift 2
while [ $# -gt 0 ]; do
	case "$1" in
		--index=*)
			INDEX_NAME="${1#*=}"
			;;
		--merge)
			MERGE_INDICES=true
			;;
		--run-mageck)
			RUN_MAGECK=true
			;;
		--mode=*)
			MAGECK_MODE="${1#*=}"
			if [ "$MAGECK_MODE" != "test" ] && [ "$MAGECK_MODE" != "mle" ]; then
				echo "错误: mageck模式必须是test或mle"
				show_usage
			fi
			;;
		*)
			echo "错误: 未知选项 $1"
			show_usage
			;;
	esac
	shift
done

BASE_DIR="/home/gyk/project/SiCRIS"
SCRIPT_DIR="${BASE_DIR}/scripts"
RESULTS_DIR="${BASE_DIR}/results"
LOG_DIR="${BASE_DIR}/logs/batch_process"
MAGECK_OUTPUT_DIR="${RESULTS_DIR}/mageck_output"
CHECKPOINT_DIR="${LOG_DIR}/checkpoints"

# 清理检查点
rm -rf "${CHECKPOINT_DIR}"/*
if [ "$RUN_MAGECK" = true ]; then
	mkdir -p "${MAGECK_OUTPUT_DIR}"
fi

# 日志记录函数
log_message() {
	local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
	echo "[${timestamp}] $1"
	echo "[${timestamp}] $1" >> "${LOG_DIR}/batch_process.log"
}

# 错误处理函数
handle_error() {
	local error_message="$1"
	local checkpoint_name="$2"
	local checkpoint_file="${CHECKPOINT_DIR}/${checkpoint_name}.done"
	log_message "错误: ${error_message}"
	log_message "删除检查点: ${checkpoint_name}"
	rm -f "${checkpoint_file}"
	log_message "请修复错误后重新运行脚本"
	exit 1
}

# 检查点函数
set_checkpoint() {
	local checkpoint_name="$1"
	local checkpoint_file="${CHECKPOINT_DIR}/${checkpoint_name}.done"
	touch "${checkpoint_file}"
	log_message "设置检查点: ${checkpoint_name}"
}

check_checkpoint() {
	local checkpoint_name="$1"
	local checkpoint_file="${CHECKPOINT_DIR}/${checkpoint_name}.done"
	if [ -f "${checkpoint_file}" ]; then
		log_message "检查点已存在: ${checkpoint_name}，跳过此步骤"
		return 0
	else
		return 1
	fi
}

# 显示脚本标题
cat << "EOF"
 ____   _   ____  ____  ___  ____     ____
/ ___| (_) / ___||  _ \|_ _|/ ___|   |  _ \ _ __ ___   ___ ___  ___ ___
\___ \ | || |    | |_) || | \___ \   | |_) | '__/ _ \ / __/ _ \/ __/ __|
 ___) || || |___ |  _ < | |  ___) |  |  __/| | | (_) | (_|  __/\__ \__ \
|____/ |_| \____||_| \_\___||____/   |_|   |_|  \___/ \___\___||___/___/
EOF

log_message "开始批量处理样本"
log_message "样本表文件: ${SAMPLESHEET}"
log_message "裁剪长度: ${CROP_LENGTH}"
if [ -n "${INDEX_NAME}" ]; then
	log_message "使用指定索引: ${INDEX_NAME}"
fi
if [ "$MERGE_INDICES" = true ]; then
	log_message "将合并gene和promoter索引数据"
fi
if [ "$RUN_MAGECK" = true ]; then
	log_message "将自动运行mageck分析，模式: ${MAGECK_MODE}"
fi

# 读取样本表并处理每个样本
log_message "读取样本表并处理每个样本"

# 获取所有样本的条件信息
declare -A sample_conditions
while IFS=, read -r sample fastq_1 fastq_2 condition; do
	# 跳过标题行和空行
	if [ -z "${sample}" ] || [ "${sample}" = "sample" ]; then
		continue
	fi
	sample_conditions["${sample}"]="${condition}"
done < "${SAMPLESHEET}"

# 处理样本
process_samples() {
	# 检查是否已完成样本处理
	if check_checkpoint "samples_processed"; then
		log_message "所有样本已处理完成，跳过样本处理步骤"
		return
	fi

	log_message "开始处理样本"

	# 跳过标题行，读取样本信息
	tail -n +2 "${SAMPLESHEET}" | while IFS=, read -r sample fastq_1 fastq_2 condition; do
		# 跳过空行
		if [ -z "${sample}" ]; then
			continue
		fi

		# 检查样本是否已处理
		if check_checkpoint "sample_${sample}"; then
			log_message "样本 ${sample} 已处理，跳过"
			continue
		fi

		log_message "处理样本: ${sample}, 条件: ${condition}"

		# 检查输入文件是否存在
		if [ ! -f "${fastq_1}" ]; then
			log_message "错误: 输入文件 ${fastq_1} 不存在，跳过样本 ${sample}"
			continue
		fi

		if [ ! -f "${fastq_2}" ]; then
			log_message "错误: 输入文件 ${fastq_2} 不存在，跳过样本 ${sample}"
			continue
		fi

		# 运行get_sgRNAcount.sh脚本处理样本
		log_message "运行get_sgRNAcount.sh脚本处理样本 ${sample}"
		if [ -n "${INDEX_NAME}" ]; then
			"${SCRIPT_DIR}/get_sgRNAcount.sh" "${fastq_1}" "${fastq_2}" "${sample}" "${CROP_LENGTH}" "${INDEX_NAME}"
		else
			"${SCRIPT_DIR}/get_sgRNAcount.sh" "${fastq_1}" "${fastq_2}" "${sample}" "${CROP_LENGTH}"
		fi

		log_message "样本 ${sample} 处理完成"
		set_checkpoint "sample_${sample}"
	done

	set_checkpoint "samples_processed"
	log_message "所有样本处理完成"
}

# 跳过标题行，读取样本信息
process_samples

# 整理结果，准备mageck输入文件
if check_checkpoint "results_prepared"; then
	log_message "结果文件已准备完成，跳过结果整理步骤"
else
	log_message "整理结果，准备mageck输入文件"

# 找出所有不同的条件
conditions=()
for condition in "${sample_conditions[@]}"; do
	if [[ ! " ${conditions[@]} " =~ " ${condition} " ]]; then
		conditions+=("${condition}")
	fi
done

log_message "检测到以下条件: ${conditions[*]}"

# 找出控制组条件
control_condition="control"
log_message "使用 '${control_condition}' 作为对照组"

# 处理基因sgRNA计数文件
if [ -n "${INDEX_NAME}" ]; then
	# 使用指定索引
	process_indices=("${INDEX_NAME}")
else
	# 使用基因和启动子索引
	process_indices=("gene" "promoter")
fi

# 处理每个索引
for process_index in "${process_indices[@]}"; do
	# 检查索引是否已处理
	if check_checkpoint "index_${process_index}_processed"; then
		log_message "索引 ${process_index} 已处理，跳过"
		continue
	fi
	log_message "处理 ${process_index} 索引的计数文件"

	# 创建临时目录
	tmp_dir="${RESULTS_DIR}/mageck_input/tmp_${process_index}"
	mkdir -p "${tmp_dir}"

	# 第一个样本的文件用于获取sgRNA和基因信息
	first_sample=""
	for sample in "${!sample_conditions[@]}"; do
		first_sample="${sample}"
		break
	done

	if [ -z "${first_sample}" ]; then
		log_message "错误: 无法找到有效的样本"
		exit 1
	fi

	# 确定计数文件路径
	if [ -n "${INDEX_NAME}" ]; then
		count_file="${RESULTS_DIR}/${first_sample}/${first_sample}.top90.count.mageck.txt"
	else
		count_file="${RESULTS_DIR}/${first_sample}/${first_sample}.${process_index}.top90.count.mageck.txt"
	fi

	# 检查计数文件是否存在
	if [ ! -f "${count_file}" ]; then
		log_message "错误: 找不到计数文件 ${count_file}"
		continue
	fi

	# 提取sgRNA和基因信息
	log_message "从样本 ${first_sample} 提取sgRNA和基因信息 (${process_index})"
	cut -f1,2 "${count_file}" | tail -n +2 > "${tmp_dir}/sgrna_gene.txt"

	# 创建合并的计数文件头部，确保样本名称在第一行
	printf "sgRNA\tGene" > "${tmp_dir}/count_header.txt"
	
	# 为每个样本添加列
	for sample in "${!sample_conditions[@]}"; do
		printf "\t%s" "${sample}" >> "${tmp_dir}/count_header.txt"
	done
	printf "\n" >> "${tmp_dir}/count_header.txt"

	# 为每个sgRNA准备一个数组来存储所有样本的计数
	declare -A sgrna_counts

	# 读取sgRNA和基因信息
	while IFS=$'\t' read -r sgrna gene; do
		sgrna_counts["${sgrna}"]=""
	done < "${tmp_dir}/sgrna_gene.txt"

	# 处理每个样本的计数
	for sample in "${!sample_conditions[@]}"; do
		log_message "处理样本 ${sample} 的计数数据 (${process_index})"

		# 确定计数文件路径
		if [ -n "${INDEX_NAME}" ]; then
			sample_count_file="${RESULTS_DIR}/${sample}/${sample}.top90.count.mageck.txt"
		else
			sample_count_file="${RESULTS_DIR}/${sample}/${sample}.${process_index}.top90.count.mageck.txt"
		fi

		# 检查计数文件是否存在
		if [ ! -f "${sample_count_file}" ]; then
			log_message "警告: 找不到样本 ${sample} 的计数文件 ${sample_count_file}"
			# 为所有sgRNA添加0计数
			for sgrna in "${!sgrna_counts[@]}"; do
				sgrna_counts["${sgrna}"]="${sgrna_counts["${sgrna}"]}\t0"
			done
			continue
		fi

		# 创建临时计数映射文件
		awk -F'\t' 'NR>1 {print $1 "\t" $3}' "${sample_count_file}" > "${tmp_dir}/${sample}_counts.txt"

		# 读取样本计数并添加到sgrna_counts
		declare -A sample_count_map
		while IFS=$'\t' read -r sgrna count; do
			sample_count_map["${sgrna}"]="${count}"
		done < "${tmp_dir}/${sample}_counts.txt"

		# 为每个sgRNA添加计数
		for sgrna in "${!sgrna_counts[@]}"; do
			if [ -n "${sample_count_map["${sgrna}"]}" ]; then
				sgrna_counts["${sgrna}"]="${sgrna_counts["${sgrna}"]}\t${sample_count_map["${sgrna}"]}"
			else
				sgrna_counts["${sgrna}"]="${sgrna_counts["${sgrna}"]}\t0"
			fi
		done
	done

	# 创建合并的计数文件
	log_message "创建合并的计数文件 (${process_index})"
	> "${tmp_dir}/counts_only.txt"

	# 将sgRNA和计数写入文件
	while IFS=$'\t' read -r sgrna gene; do
		echo -e "${sgrna}\t${gene}${sgrna_counts["${sgrna}"]}" >> "${tmp_dir}/counts_only.txt"
	done < "${tmp_dir}/sgrna_gene.txt"

	# 合并头部和计数
	cat "${tmp_dir}/count_header.txt" "${tmp_dir}/counts_only.txt" > "${RESULTS_DIR}/mageck_input/all_samples_${process_index}.count.txt"

	log_message "创建了合并的计数文件: ${RESULTS_DIR}/mageck_input/all_samples_${process_index}.count.txt"
	set_checkpoint "index_${process_index}_processed"

	# 为每个非控制条件创建样本索引文件
	for condition in "${conditions[@]}"; do
		if [ "${condition}" != "${control_condition}" ]; then
			log_message "为条件 '${condition}' 创建样本索引文件 (${process_index})"

			# 创建样本索引文件
			echo -e "sample\tcontrol\ttreatment" > "${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_${process_index}.txt"

			# 添加样本信息
			for sample in "${!sample_conditions[@]}"; do
				if [ "${sample_conditions[${sample}]}" = "${control_condition}" ]; then
					echo -e "${sample}\t1\t0" >> "${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_${process_index}.txt"
				elif [ "${sample_conditions[${sample}]}" = "${condition}" ]; then
					echo -e "${sample}\t0\t1" >> "${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_${process_index}.txt"
				fi
			done

			log_message "创建了样本索引文件: ${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_${process_index}.txt"
		fi
	done
done

# 如果需要合并gene和promoter数据
if [ "$MERGE_INDICES" = true ] && [ -z "${INDEX_NAME}" ]; then
	# 检查是否已合并数据
	if check_checkpoint "indices_merged"; then
		log_message "gene和promoter数据已合并，跳过合并步骤"
	else
	log_message "合并gene和promoter索引数据"

	# 创建临时目录
	tmp_dir="${RESULTS_DIR}/mageck_input/tmp_merged"
	mkdir -p "${tmp_dir}"

	# 检查gene和promoter计数文件是否存在
	gene_count_file="${RESULTS_DIR}/mageck_input/all_samples_gene.count.txt"
	promoter_count_file="${RESULTS_DIR}/mageck_input/all_samples_promoter.count.txt"

	if [ ! -f "${gene_count_file}" ] || [ ! -f "${promoter_count_file}" ]; then
		log_message "错误: 找不到gene或promoter计数文件，无法合并"
	else
		# 合并gene和promoter计数文件
		log_message "合并gene和promoter计数文件"

		# 检查并修复gene计数文件的格式
		log_message "检查并修复gene计数文件的格式"
		# 提取第一行作为列标题，第二行作为样本名称，合并成正确的格式
		awk 'BEGIN {FS="\t"; OFS="\t"} 
		     NR==1 {header=$0} 
		     NR==2 {samples=$0; split(samples,arr,"\t"); for(i=1;i<=NF;i++) {if(i==1) {printf "%s\t%s", $1, $2} else if(i>2) {printf "\t%s", arr[i-2]}}; printf "\n"} 
		     NR>2 {print $0}' "${gene_count_file}" > "${tmp_dir}/gene_count_fixed.txt"

		# 检查并修复promoter计数文件的格式
		log_message "检查并修复promoter计数文件的格式"
		# 提取第一行作为列标题，第二行作为样本名称，合并成正确的格式
		awk 'BEGIN {FS="\t"; OFS="\t"} 
		     NR==1 {header=$0} 
		     NR==2 {samples=$0; split(samples,arr,"\t"); for(i=1;i<=NF;i++) {if(i==1) {printf "%s\t%s", $1, $2} else if(i>2) {printf "\t%s", arr[i-2]}}; printf "\n"} 
		     NR>2 {print $0}' "${promoter_count_file}" > "${tmp_dir}/promoter_count_fixed.txt"

		# 提取头部行
		head -n 1 "${tmp_dir}/gene_count_fixed.txt" > "${tmp_dir}/merged_header.txt"

		# 合并数据行（跳过头部）
		tail -n +2 "${tmp_dir}/gene_count_fixed.txt" > "${tmp_dir}/gene_data.txt"
		tail -n +2 "${tmp_dir}/promoter_count_fixed.txt" > "${tmp_dir}/promoter_data.txt"

		# 合并数据
		cat "${tmp_dir}/gene_data.txt" "${tmp_dir}/promoter_data.txt" > "${tmp_dir}/merged_data.txt"

		# 创建合并的计数文件
		cat "${tmp_dir}/merged_header.txt" "${tmp_dir}/merged_data.txt" > "${RESULTS_DIR}/mageck_input/all_samples_merged.count.txt"

		# 检查合并后的文件格式是否正确
		log_message "检查合并后的文件格式"
		# 显示文件头部以进行调试
		log_message "文件头部: $(head -n 1 "${RESULTS_DIR}/mageck_input/all_samples_merged.count.txt")"
		# 修复文件格式，删除可能的转义字符
		sed -i '1s/\\t/\t/g' "${RESULTS_DIR}/mageck_input/all_samples_merged.count.txt"
		# 再次检查格式
		if ! head -n 1 "${RESULTS_DIR}/mageck_input/all_samples_merged.count.txt" | grep -q "piRNA-Rep1"; then
			log_message "尝试使用备用方法修复文件格式"
			# 备用方法：直接创建正确的头部
			printf "sgRNA\tGene" > "${tmp_dir}/fixed_header.txt"
			for sample in "${!sample_conditions[@]}"; do
				printf "\t%s" "${sample}" >> "${tmp_dir}/fixed_header.txt"
			done
			printf "\n" >> "${tmp_dir}/fixed_header.txt"
			# 将新头部与原数据合并
			tail -n +2 "${RESULTS_DIR}/mageck_input/all_samples_merged.count.txt" > "${tmp_dir}/merged_data_only.txt"
			cat "${tmp_dir}/fixed_header.txt" "${tmp_dir}/merged_data_only.txt" > "${RESULTS_DIR}/mageck_input/all_samples_merged.count.txt"
			# 再次检查
			if ! head -n 1 "${RESULTS_DIR}/mageck_input/all_samples_merged.count.txt" | grep -q "piRNA-Rep1"; then
				handle_error "合并后的文件格式不正确，样本名称未正确显示在标题行" "indices_merged"
			fi
		fi

		log_message "创建了合并的计数文件: ${RESULTS_DIR}/mageck_input/all_samples_merged.count.txt"
		set_checkpoint "indices_merged"

		# 为每个非控制条件创建样本索引文件
		for condition in "${conditions[@]}"; do
			if [ "${condition}" != "${control_condition}" ]; then
				log_message "为条件 '${condition}' 创建样本索引文件 (merged)"

				# 创建样本索引文件
				echo -e "sample\tcontrol\ttreatment" > "${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_merged.txt"

				# 添加样本信息
				for sample in "${!sample_conditions[@]}"; do
					if [ "${sample_conditions[${sample}]}" = "${control_condition}" ]; then
						echo -e "${sample}\t1\t0" >> "${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_merged.txt"
					elif [ "${sample_conditions[${sample}]}" = "${condition}" ]; then
						echo -e "${sample}\t0\t1" >> "${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_merged.txt"
					fi
				done

				log_message "创建了样本索引文件: ${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_merged.txt"
				set_checkpoint "design_${condition}_vs_${control_condition}_merged"
			fi
		done
		fi
	fi
fi

# 清理临时文件
if ! check_checkpoint "tmp_files_cleaned"; then
	log_message "清理临时文件"
	rm -rf "${RESULTS_DIR}/mageck_input/tmp"*
	set_checkpoint "tmp_files_cleaned"
fi

set_checkpoint "results_prepared"

log_message "批量处理完成"
echo "批量处理完成，结果文件位于: ${RESULTS_DIR}/mageck_input/"

# 如果需要自动运行mageck
if [ "$RUN_MAGECK" = true ]; then
	if check_checkpoint "mageck_analysis_completed"; then
		log_message "mageck分析已完成，跳过分析步骤"
	else
		log_message "开始运行mageck分析，模式: ${MAGECK_MODE}"

	# 确定要处理的索引
	if [ -n "${INDEX_NAME}" ]; then
		run_indices=("${INDEX_NAME}")
	elif [ "$MERGE_INDICES" = true ]; then
		run_indices=("merged")
	else
		run_indices=("gene" "promoter")
	fi

	# 为每个索引和条件运行mageck
		for run_index in "${run_indices[@]}"; do
			for condition in "${conditions[@]}"; do
				if [ "${condition}" != "${control_condition}" ]; then
					# 检查是否已完成此条件的mageck分析
					if check_checkpoint "mageck_${MAGECK_MODE}_${condition}_vs_${control_condition}_${run_index}"; then
						log_message "mageck ${MAGECK_MODE} 分析已完成: ${condition} vs ${control_condition} (${run_index})，跳过"
						continue
					fi
				log_message "运行mageck ${MAGECK_MODE} 分析: ${condition} vs ${control_condition} (${run_index})"

				# 确定输入文件
				count_file="${RESULTS_DIR}/mageck_input/all_samples_${run_index}.count.txt"
				design_file="${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_${run_index}.txt"
				output_prefix="${MAGECK_OUTPUT_DIR}/${condition}_vs_${control_condition}_${run_index}"

				# 检查输入文件是否存在
				if [ ! -f "${count_file}" ] || [ ! -f "${design_file}" ]; then
					log_message "错误: 找不到输入文件，跳过mageck分析"
					continue
				fi

				# 运行mageck
					if [ "${MAGECK_MODE}" = "test" ]; then
						log_message "运行mageck test: ${output_prefix}"
						mageck test -k "${count_file}" -d "${design_file}" -n "${output_prefix}" --pdf-report || handle_error "mageck test分析失败" "mageck_${MAGECK_MODE}_${condition}_vs_${control_condition}_${run_index}"
					else
						log_message "运行mageck mle: ${output_prefix}"
						mageck mle -k "${count_file}" -d "${design_file}" -n "${output_prefix}" || handle_error "mageck mle分析失败" "mageck_${MAGECK_MODE}_${condition}_vs_${control_condition}_${run_index}"
					fi

				log_message "mageck ${MAGECK_MODE} 分析完成: ${output_prefix}"
					set_checkpoint "mageck_${MAGECK_MODE}_${condition}_vs_${control_condition}_${run_index}"
			fi
		done
	done

	log_message "所有mageck分析完成"
		set_checkpoint "mageck_analysis_completed"
		echo "所有mageck分析完成，结果文件位于: ${MAGECK_OUTPUT_DIR}/"
	fi
else
	echo "可以使用以下命令运行mageck分析:"
	echo ""

	# 显示mageck命令示例
	if [ -n "${INDEX_NAME}" ]; then
		for condition in "${conditions[@]}"; do
			if [ "${condition}" != "${control_condition}" ]; then
				echo "mageck test -k ${RESULTS_DIR}/mageck_input/all_samples_${INDEX_NAME}.count.txt -d ${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_${INDEX_NAME}.txt -n ${MAGECK_OUTPUT_DIR}/${condition}_vs_${control_condition}_${INDEX_NAME}"
				echo "mageck mle -k ${RESULTS_DIR}/mageck_input/all_samples_${INDEX_NAME}.count.txt -d ${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_${INDEX_NAME}.txt -n ${MAGECK_OUTPUT_DIR}/${condition}_vs_${control_condition}_${INDEX_NAME}"
			fi
		done
	elif [ "$MERGE_INDICES" = true ]; then
		for condition in "${conditions[@]}"; do
			if [ "${condition}" != "${control_condition}" ]; then
				echo "mageck test -k ${RESULTS_DIR}/mageck_input/all_samples_merged.count.txt -d ${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_merged.txt -n ${MAGECK_OUTPUT_DIR}/${condition}_vs_${control_condition}_merged"
				echo "mageck mle -k ${RESULTS_DIR}/mageck_input/all_samples_merged.count.txt -d ${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_merged.txt -n ${MAGECK_OUTPUT_DIR}/${condition}_vs_${control_condition}_merged"
			fi
		done
	else
		for process_index in "gene" "promoter"; do
			for condition in "${conditions[@]}"; do
				if [ "${condition}" != "${control_condition}" ]; then
					echo "mageck test -k ${RESULTS_DIR}/mageck_input/all_samples_${process_index}.count.txt -d ${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_${process_index}.txt -n ${MAGECK_OUTPUT_DIR}/${condition}_vs_${control_condition}_${process_index}"
					echo "mageck mle -k ${RESULTS_DIR}/mageck_input/all_samples_${process_index}.count.txt -d ${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_${process_index}.txt -n ${MAGECK_OUTPUT_DIR}/${condition}_vs_${control_condition}_${process_index}"
				fi
			done
		done
		fi
	fi
fi

exit 0