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
	echo "  --run-mageck   - 自动运行mageck分析"
	echo "  --mode=MODE    - mageck分析模式 (test或mle，默认为test)"
	echo "示例: $0 /home/gyk/project/SiCRIS/data/samplesheet.csv 100 --run-mageck --mode=mle"
	exit 1
}

# 参数检查
if [ $# -lt 2 ]; then
	echo "错误: 参数不足"
	show_usage
fi

# 检查必要的工具
if ! command -v csvtk &> /dev/null; then
	echo "错误: 未找到 csvtk 工具，请先安装 csvtk"
	echo "安装方法: conda install -c bioconda csvtk 或 pip install csvtk"
	exit 1
fi

# 检查样本表文件是否存在
if [ ! -f "$1" ]; then
	echo "错误: 样本表文件 $1 不存在"
	exit 1
fi

SAMPLESHEET="$1"
CROP_LENGTH="$2"
INDEX_NAME=""
RUN_MAGECK=false
MAGECK_MODE="test"

# 解析可选参数
shift 2
while [ $# -gt 0 ]; do
	case "$1" in
		--index=*)
			INDEX_NAME="${1#*=}"
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

# 处理基因和启动子sgRNA计数文件
process_indices=("gene" "promoter")

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
	
	# 使用csvtk提取sgRNA和基因信息
	csvtk cut -t -f 1,2 "${count_file}" | csvtk del-header -t > "${tmp_dir}/sgrna_gene.txt"
	
	# 创建所有样本的计数文件列表
	sample_count_files=()
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
			continue
		fi
		
		# 提取sgRNA和计数信息，并重命名列
		csvtk cut -t -f 1,3 "${sample_count_file}" | csvtk del-header -t | \
		csvtk add-header -t -n "sgRNA,${sample}" > "${tmp_dir}/${sample}_counts.txt"
		
		sample_count_files+=("${tmp_dir}/${sample}_counts.txt")
	done
	
	# 使用csvtk合并所有样本的计数文件
	log_message "合并所有样本的计数数据 (${process_index})"
	
	# 首先创建带有sgRNA和基因信息的基础文件
	csvtk add-header -t -n "sgRNA,Gene" "${tmp_dir}/sgrna_gene.txt" > "${tmp_dir}/base_file.txt"
	
	# 逐个合并样本计数文件
	merged_file="${tmp_dir}/base_file.txt"
	for sample_file in "${sample_count_files[@]}"; do
		temp_merged="${tmp_dir}/temp_merged_$RANDOM.txt"
		csvtk join -t -f "sgRNA" "$merged_file" "$sample_file" > "$temp_merged"
		merged_file="$temp_merged"
	done
	
	# 处理缺失值，将NA替换为0
	csvtk replace -t -p "NA" -r "0" "$merged_file" > "${RESULTS_DIR}/mageck_input/all_samples_${process_index}.count.txt"

	log_message "创建了合并的计数文件: ${RESULTS_DIR}/mageck_input/all_samples_${process_index}.count.txt"
	set_checkpoint "index_${process_index}_processed"

	# 为每个非控制条件创建样本索引文件
	for condition in "${conditions[@]}"; do
		if [ "${condition}" != "${control_condition}" ]; then
			log_message "为条件 '${condition}' 创建样本索引文件 (${process_index})"

			# 创建样本索引文件
			design_file="${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_${process_index}.txt"
			
			# 使用csvtk创建设计文件
			echo -e "sample\tcontrol\ttreatment" > "$design_file"
			
			# 使用csvtk添加样本信息
			for sample in "${!sample_conditions[@]}"; do
				if [ "${sample_conditions[${sample}]}" = "${control_condition}" ]; then
					echo -e "${sample}\t1\t0" >> "$design_file"
				elif [ "${sample_conditions[${sample}]}" = "${condition}" ]; then
					echo -e "${sample}\t0\t1" >> "$design_file"
				fi
			done

			log_message "创建了样本索引文件: ${RESULTS_DIR}/mageck_input/${condition}_vs_${control_condition}_${process_index}.txt"
		fi
	done
done



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