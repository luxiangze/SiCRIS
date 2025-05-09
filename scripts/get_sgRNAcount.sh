#!/bin/bash
#
# SiCRIS sgRNA计数脚本
# 用途：处理sgRNA测序数据，进行比对和计数
#
# 用法: $0 <R1.fq.gz> <R2.fq.gz> <sample_name> <crop_length> [index_name]
#   R1.fq.gz      - 正向测序数据文件
#   R2.fq.gz      - 反向测序数据文件
#   sample_name   - 样本名称
#   crop_length   - 裁剪长度
#   index_name    - 可选，索引名称
#

set -e  # 遇到错误立即退出
set -o pipefail  # 管道中任何命令失败都会导致整个管道失败

# 显示使用方法
function show_usage {
    echo "用法: $0 <R1.fq.gz> <R2.fq.gz> <sample_name> <crop_length> [index_name]"
    echo "示例: $0 /path/to/sample_R1.fq.gz /path/to/sample_R2.fq.gz sample1 100"
    exit 1
}

# 参数检查
if [ $# -lt 4 ]; then
    echo "错误: 参数不足"
    show_usage
fi

# 检查输入文件是否存在
if [ ! -f "$1" ]; then
    echo "错误: 输入文件 $1 不存在"
    exit 1
fi

if [ ! -f "$2" ]; then
    echo "错误: 输入文件 $2 不存在"
    exit 1
fi

# 定义基础路径
BASE_DIR="/home/gyk/project/SiCRIS"
INDEX_DIR="/home/gyk/project/SiCRIS/data/index"

# 设置目录和路径变量
SAMPLE_NAME="$3"
LOG_DIR="${BASE_DIR}/logs/${SAMPLE_NAME}"
WORK_DIR="${BASE_DIR}/work/${SAMPLE_NAME}"
RESULT_DIR="${BASE_DIR}/results/${SAMPLE_NAME}"

# 创建必要的目录
mkdir -p "${LOG_DIR}"
mkdir -p "${WORK_DIR}"
mkdir -p "${RESULT_DIR}"

# 日志记录函数
log_message() {
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[${timestamp}] $1"
    echo "[${timestamp}] $1" >> "${LOG_DIR}/${SAMPLE_NAME}_main.log"
}

log_section_start() {
    echo "**************** $1 ****************"
    log_message "开始: $1"
    echo "*************************************************************"
}

log_section_end() {
    echo "**************** $1 ****************"
    log_message "完成: $1"
    echo "*************************************************************"
    # 创建检查点标记文件
    touch "${WORK_DIR}/.${SAMPLE_NAME}_$(echo $1 | tr ' ' '_')_done"
}

# 检查步骤是否已完成
check_step_done() {
    local step_name=$(echo $1 | tr ' ' '_')
    if [ -f "${WORK_DIR}/.${SAMPLE_NAME}_${step_name}_done" ]; then
        log_message "检测到步骤 '$1' 已完成，跳过此步骤"
        return 0  # 步骤已完成
    else
        return 1  # 步骤未完成
    fi
}

# 显示脚本标题
cat << "EOF"
 ____   _   ____  ____  ___  ____
/ ___| (_) / ___||  _ \|_ _|/ ___|
\___ \ | || |    | |_) || | \___ \
 ___) || || |___ |  _ < | |  ___) |
|____/ |_| \____||_| \_\___||____/
EOF

log_message "开始处理样本: ${SAMPLE_NAME}"

# 步骤1: 修剪PE读数
if check_step_done "修剪PE读数"; then
    # 如果步骤已完成，检查所需的文件是否存在
    if [ ! -f "${WORK_DIR}/${SAMPLE_NAME}_clean_R1.fq.gz" ] || [ ! -f "${WORK_DIR}/${SAMPLE_NAME}_clean_R2.fq.gz" ]; then
        log_message "警告: 修剪后的文件不存在，将重新执行修剪步骤"
        rm -f "${WORK_DIR}/.${SAMPLE_NAME}_修剪PE读数_done"
    else
        log_message "修剪后的文件已存在，跳过修剪步骤"
    fi
fi

# 如果步骤未完成，执行修剪操作
if ! check_step_done "修剪PE读数"; then
    log_section_start "修剪PE读数"

    # 创建输入文件的软链接
    ln -sf "$1" "${WORK_DIR}/${SAMPLE_NAME}_R1.fq.gz"
    ln -sf "$2" "${WORK_DIR}/${SAMPLE_NAME}_R2.fq.gz"

    # 使用trimmomatic进行质量修剪
    log_message "运行trimmomatic进行质量修剪，裁剪长度: $4"
    trimmomatic PE -threads 20 \
        "${WORK_DIR}/${SAMPLE_NAME}_R1.fq.gz" \
        "${WORK_DIR}/${SAMPLE_NAME}_R2.fq.gz" \
        "${WORK_DIR}/${SAMPLE_NAME}_clean_R1.fq.gz" \
        "${WORK_DIR}/${SAMPLE_NAME}_un_R1.fq.gz" \
        "${WORK_DIR}/${SAMPLE_NAME}_clean_R2.fq.gz" \
        "${WORK_DIR}/${SAMPLE_NAME}_un_R2.fq.gz" \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 CROP:$4 \
        2>"${LOG_DIR}/${SAMPLE_NAME}_01_trim.log"

    # 删除未配对的读数文件
    rm -f "${WORK_DIR}"/*_un_*

    log_section_end "修剪PE读数"
fi
# 步骤2: 使用Flash合并PE读数
if check_step_done "Flash合并PE读数"; then
    # 如果步骤已完成，检查所需的文件是否存在
    if [ ! -f "${WORK_DIR}/${SAMPLE_NAME}.extendedFrags.cutadapt.fastq.gz" ]; then
        log_message "警告: 合并后的文件不存在，将重新执行合并步骤"
        rm -f "${WORK_DIR}/.${SAMPLE_NAME}_Flash合并PE读数_done"
    else
        log_message "合并后的文件已存在，跳过合并步骤"
    fi
fi

# 如果步骤未完成，执行合并操作
if ! check_step_done "Flash合并PE读数"; then
    log_section_start "Flash合并PE读数"

    # 使用Flash合并配对的读数
    log_message "运行Flash合并配对的读数"
    # 确保工作目录存在
    mkdir -p "${WORK_DIR}"

    # 检查修剪后的文件是否存在
    if [ ! -f "${WORK_DIR}/${SAMPLE_NAME}_clean_R1.fq.gz" ] || [ ! -f "${WORK_DIR}/${SAMPLE_NAME}_clean_R2.fq.gz" ]; then
        log_message "错误: 修剪后的文件不存在，无法执行合并步骤"
        exit 1
    fi

    # 切换到工作目录运行flash，避免路径问题
    cd "${WORK_DIR}" && \
    flash -z -O -M 149 \
        "${SAMPLE_NAME}_clean_R1.fq.gz" \
        "${SAMPLE_NAME}_clean_R2.fq.gz" \
        -o "${SAMPLE_NAME}" \
        > "${LOG_DIR}/${SAMPLE_NAME}_02_flash.log" || {
            log_message "错误: Flash处理失败"
            exit 1
        }

    # 切回原目录
    cd - > /dev/null

    # 删除不需要的中间文件
    rm -f "${WORK_DIR}"/*.hist* "${WORK_DIR}"/*.notCombined*

    # 使用cutadapt进行序列适配器剔除
    log_message "运行cutadapt进行序列适配器剔除"

    # 切换到工作目录运行cutadapt，避免路径问题
    cd "${WORK_DIR}" && \
    cutadapt -g TAGCTCTAAAAC -g GCTCTACAAGTG \
        --length 20 \
        --discard-untrimmed \
        -m 20 \
        --cores=20 \
        -o "${SAMPLE_NAME}.extendedFrags.cutadapt.fastq.gz" \
        "${SAMPLE_NAME}.extendedFrags.fastq.gz" \
        > "${LOG_DIR}/${SAMPLE_NAME}_03_cutadapt.log" || {
            log_message "错误: Cutadapt处理失败"
            exit 1
        }

    # 切回原目录
    cd - > /dev/null

    log_section_end "Flash合并PE读数"
fi

# 步骤3: 使用Bowtie2进行sgRNA比对
if check_step_done "Bowtie2 sgRNA比对"; then
    # 如果步骤已完成，检查所需的文件是否存在
    if [ -n "$5" ]; then
        # 指定索引模式
        if [ ! -f "${WORK_DIR}/${SAMPLE_NAME}.extendedFrags.cutadapt.fastq.$5.sam" ]; then
            log_message "警告: 指定索引的比对文件不存在，将重新执行比对步骤"
            rm -f "${WORK_DIR}/.${SAMPLE_NAME}_Bowtie2_sgRNA比对_done"
        else
            log_message "指定索引的比对文件已存在，跳过比对步骤"
        fi
    else
        # 默认模式（基因和启动子）
        if [ ! -f "${WORK_DIR}/${SAMPLE_NAME}.extendedFrags.cutadapt.fastq.gene.sam" ] ||
           [ ! -f "${WORK_DIR}/${SAMPLE_NAME}.extendedFrags.cutadapt.fastq.promoter.sam" ]; then
            log_message "警告: 基因或启动子比对文件不存在，将重新执行比对步骤"
            rm -f "${WORK_DIR}/.${SAMPLE_NAME}_Bowtie2_sgRNA比对_done"
        else
            log_message "基因和启动子比对文件已存在，跳过比对步骤"
        fi
    fi
fi

# 如果步骤未完成，执行比对操作
if ! check_step_done "Bowtie2 sgRNA比对"; then
    log_section_start "Bowtie2 sgRNA比对"

    # 设置公共的bowtie2参数
    BOWTIE2_PARAMS="--no-unal --no-head -p 20"
    INPUT_FASTQ="${WORK_DIR}/${SAMPLE_NAME}.extendedFrags.cutadapt.fastq.gz"

    # 检查输入文件是否存在
    if [ ! -f "${INPUT_FASTQ}" ]; then
        log_message "错误: 输入文件 ${INPUT_FASTQ} 不存在，无法执行比对步骤"
        exit 1
    fi

    # 检查是否提供了索引名称
    if [ -n "$5" ]; then
        # 使用指定的索引
        INDEX_NAME="$5"
        log_message "使用指定的索引: ${INDEX_NAME}"

        # 检查索引文件是否存在
        if [ ! -f "${INDEX_DIR}/${INDEX_NAME}.1.bt2" ] && [ ! -f "${INDEX_DIR}/${INDEX_NAME}.1.bt2l" ]; then
            log_message "错误: 找不到索引文件 ${INDEX_DIR}/${INDEX_NAME}"
            exit 1
        fi

        # 运行bowtie2比对
        log_message "运行bowtie2比对到索引: ${INDEX_NAME}"
        bowtie2 ${BOWTIE2_PARAMS} \
            -x "${INDEX_DIR}/${INDEX_NAME}" \
            -U "${INPUT_FASTQ}" \
            -S "${WORK_DIR}/${SAMPLE_NAME}.extendedFrags.cutadapt.fastq.${INDEX_NAME}.sam" \
            2>"${LOG_DIR}/${SAMPLE_NAME}_04_bowtie2.${INDEX_NAME}.log" || {
                log_message "错误: Bowtie2比对失败"
                exit 1
            }

        # 记录比对结果
        GENE_SAM_FILE="${WORK_DIR}/${SAMPLE_NAME}.extendedFrags.cutadapt.fastq.${INDEX_NAME}.sam"
        log_message "生成比对文件: ${GENE_SAM_FILE}"

    else
        # 使用默认的基因和启动子索引
        GENE_INDEX="20200725.gene.6.final"
        PROMOTER_INDEX="20200726.Promoter.3.final"

        log_message "使用默认的基因索引: ${GENE_INDEX}"
        # 检查基因索引文件是否存在
        if [ ! -f "${INDEX_DIR}/${GENE_INDEX}.1.bt2" ] && [ ! -f "${INDEX_DIR}/${GENE_INDEX}.1.bt2l" ]; then
            log_message "错误: 找不到基因索引文件 ${INDEX_DIR}/${GENE_INDEX}"
            exit 1
        fi

        # 运行bowtie2比对到基因索引
        bowtie2 ${BOWTIE2_PARAMS} \
            -x "${INDEX_DIR}/${GENE_INDEX}" \
            -U "${INPUT_FASTQ}" \
            -S "${WORK_DIR}/${SAMPLE_NAME}.extendedFrags.cutadapt.fastq.gene.sam" \
            2>"${LOG_DIR}/${SAMPLE_NAME}_04_bowtie2.gene.log" || {
                log_message "错误: 基因索引Bowtie2比对失败"
                exit 1
            }

        log_message "使用默认的启动子索引: ${PROMOTER_INDEX}"
        # 检查启动子索引文件是否存在
        if [ ! -f "${INDEX_DIR}/${PROMOTER_INDEX}.1.bt2" ] && [ ! -f "${INDEX_DIR}/${PROMOTER_INDEX}.1.bt2l" ]; then
            log_message "错误: 找不到启动子索引文件 ${INDEX_DIR}/${PROMOTER_INDEX}"
            exit 1
        fi

        # 运行bowtie2比对到启动子索引
        bowtie2 ${BOWTIE2_PARAMS} \
            -x "${INDEX_DIR}/${PROMOTER_INDEX}" \
            -U "${INPUT_FASTQ}" \
            -S "${WORK_DIR}/${SAMPLE_NAME}.extendedFrags.cutadapt.fastq.promoter.sam" \
            2>"${LOG_DIR}/${SAMPLE_NAME}_04_bowtie2.promoter.log" || {
                log_message "错误: 启动子索引Bowtie2比对失败"
                exit 1
            }

        # 记录比对结果
        log_message "生成基因和启动子比对文件"
    fi

    log_section_end "Bowtie2 sgRNA比对"
fi


# 步骤4: sgRNA计数
# 定义处理函数
process_sam_file() {
    local sam_file=$1
    local output_prefix=$2
    local index_name=$3
    local temp_prefix="${WORK_DIR}/temp_${SAMPLE_NAME}"

    log_message "处理SAM文件: ${sam_file}"

    # 提取匹配的sgRNA序列
    cat "${sam_file}" | \
        awk '{print $3,$6}' | \
        egrep '19M|20M' | \
        awk '{print $1}' | \
        sort | uniq -c | \
        awk '{print $2"\t"$1}' | \
        sort -k2nr > "${temp_prefix}_1"

    # 计算前90%的记录数
    wc -l "${temp_prefix}_1" | \
        awk '{print $1*0.9}' | \
        sed 's/\..*//g' > "${temp_prefix}_2"

    # 提取前90%的记录
    head -$(cat "${temp_prefix}_2") "${temp_prefix}_1" > "${temp_prefix}_3"

    # 生成top90计数文件
    cat "${INDEX_DIR}/${index_name}.name" "${temp_prefix}_3" | \
        sort | \
        awk '{sum[$1]+=$2}END{for(c in sum){print c,sum[c]}}' | \
        sed 's/_BmEg/\tBmEg/g' | \
        sort -k1 | \
        sed 's/ /\t/g' > "${WORK_DIR}/${output_prefix}.top90.count.mageck.txt"

    # 生成原始计数文件
    cat "${INDEX_DIR}/${index_name}.name" "${temp_prefix}_1" | \
        sort | \
        awk '{sum[$1]+=$2}END{for(c in sum){print c,sum[c]}}' | \
        sed 's/_BmEg/\tBmEg/g' | \
        sort -k1 | \
        sed 's/ /\t/g' > "${WORK_DIR}/${output_prefix}.raw.count.txt"

    # 添加标题行
    sed -i "1s/^/sgRNA\tGene\t${SAMPLE_NAME}\n/" "${WORK_DIR}/${output_prefix}.top90.count.mageck.txt"
    sed -i "1s/^/sgRNA\tGene\t${SAMPLE_NAME}\n/" "${WORK_DIR}/${output_prefix}.raw.count.txt"

    # 生成统计信息
    less "${WORK_DIR}/${output_prefix}.raw.count.txt" | \
        grep -w -v 0 | \
        sed 1d | \
        awk '{print $2}' | \
        sort | uniq -c | \
        awk '{print $1}' | \
        sort | uniq -c | \
        awk '{print $2"\t"$1}' | \
        sed "1s/^/Number_of_sgRNAs_per_gene\tNumber_of_genes\n/" > "${WORK_DIR}/${output_prefix}.sgRNAnumber.stats.txt"

    # 添加总测序基因数
    less "${WORK_DIR}/${output_prefix}.raw.count.txt" | \
        grep -w -v 0 | \
        sed 1d | \
        awk '{print $2}' | \
        sort | uniq -c | \
        wc -l | \
        awk '{print "total_sequenced_gene_number\t"$0}' >> "${WORK_DIR}/${output_prefix}.sgRNAnumber.stats.txt"

    # 添加总设计基因数
    less "${WORK_DIR}/${output_prefix}.raw.count.txt" | \
        sed 1d | \
        awk '{print $2}' | \
        sort | uniq -c | \
        wc -l | \
        awk '{print "total_desigened_gene_number\t"$0}' >> "${WORK_DIR}/${output_prefix}.sgRNAnumber.stats.txt"

    # 删除临时文件
    rm -f "${temp_prefix}_"*

    log_message "生成统计文件: ${WORK_DIR}/${output_prefix}.sgRNAnumber.stats.txt"
}

# 检查是否已完成sgRNA计数步骤
if check_step_done "sgRNA计数"; then
    # 如果步骤已完成，检查所需的文件是否存在
    if [ -n "$5" ]; then
        # 指定索引模式
        if [ ! -f "${WORK_DIR}/${SAMPLE_NAME}.raw.count.txt" ] ||
           [ ! -f "${WORK_DIR}/${SAMPLE_NAME}.sgRNAnumber.stats.txt" ]; then
            log_message "警告: 计数结果文件不存在，将重新执行计数步骤"
            rm -f "${WORK_DIR}/.${SAMPLE_NAME}_sgRNA计数_done"
        else
            log_message "计数结果文件已存在，跳过计数步骤"
        fi
    else
        # 默认模式（基因和启动子）
        if [ ! -f "${WORK_DIR}/${SAMPLE_NAME}.gene.raw.count.txt" ] ||
           [ ! -f "${WORK_DIR}/${SAMPLE_NAME}.promoter.raw.count.txt" ] ||
           [ ! -f "${WORK_DIR}/${SAMPLE_NAME}.gene.sgRNAnumber.stats.txt" ] ||
           [ ! -f "${WORK_DIR}/${SAMPLE_NAME}.promoter.sgRNAnumber.stats.txt" ]; then
            log_message "警告: 计数结果文件不存在，将重新执行计数步骤"
            rm -f "${WORK_DIR}/.${SAMPLE_NAME}_sgRNA计数_done"
        else
            log_message "计数结果文件已存在，跳过计数步骤"
        fi
    fi
fi

# 如果步骤未完成，执行计数操作
if ! check_step_done "sgRNA计数"; then
    log_section_start "sgRNA计数"

    # 检查是否提供了索引名称
    if [ -n "$5" ]; then
        # 处理指定索引的SAM文件
        SAM_FILE="${WORK_DIR}/${SAMPLE_NAME}.extendedFrags.cutadapt.fastq.$5.sam"
        log_message "处理指定索引的SAM文件: ${SAM_FILE}"

        # 检查SAM文件是否存在
        if [ ! -f "${SAM_FILE}" ]; then
            log_message "错误: SAM文件不存在 ${SAM_FILE}"
            exit 1
        fi

        # 处理SAM文件
        process_sam_file "${SAM_FILE}" "${SAMPLE_NAME}" "$5"

        # 删除SAM文件
        rm -f "${SAM_FILE}"
        log_message "删除SAM文件以节省空间"


else
    # 处理基因和启动子索引
    GENE_SAM_FILE="${WORK_DIR}/${SAMPLE_NAME}.extendedFrags.cutadapt.fastq.gene.sam"
    PROMOTER_SAM_FILE="${WORK_DIR}/${SAMPLE_NAME}.extendedFrags.cutadapt.fastq.promoter.sam"
    GENE_INDEX_NAME="20200725.gene.6.final"
    PROMOTER_INDEX_NAME="20200726.Promoter.3.final"

    # 检查基因SAM文件是否存在
    if [ ! -f "${GENE_SAM_FILE}" ]; then
        log_message "错误: 基因SAM文件不存在 ${GENE_SAM_FILE}"
        exit 1
    fi

    # 检查启动子SAM文件是否存在
    if [ ! -f "${PROMOTER_SAM_FILE}" ]; then
        log_message "错误: 启动子SAM文件不存在 ${PROMOTER_SAM_FILE}"
        exit 1
    fi

    # 检查索引文件是否存在
    if [ ! -f "${INDEX_DIR}/${GENE_INDEX_NAME}.name" ]; then
        log_message "错误: 基因索引名称文件不存在 ${INDEX_DIR}/${GENE_INDEX_NAME}.name"
        exit 1
    fi

    if [ ! -f "${INDEX_DIR}/${PROMOTER_INDEX_NAME}.name" ]; then
        log_message "错误: 启动子索引名称文件不存在 ${INDEX_DIR}/${PROMOTER_INDEX_NAME}.name"
        exit 1
    fi

    # 处理基因SAM文件
    log_message "处理基因SAM文件"
    process_sam_file "${GENE_SAM_FILE}" "${SAMPLE_NAME}.gene" "${GENE_INDEX_NAME}"

    # 处理启动子SAM文件
    log_message "处理启动子SAM文件"
    process_sam_file "${PROMOTER_SAM_FILE}" "${SAMPLE_NAME}.promoter" "${PROMOTER_INDEX_NAME}"

    # # 删除SAM文件
    # rm -f "${GENE_SAM_FILE}" "${PROMOTER_SAM_FILE}"
    # log_message "删除SAM文件以节省空间"
	fi

    log_section_end "sgRNA计数"
fi

# 注释掉的R绘图部分，保留以便将来可能的使用
# log_section_start "R绘图"
# if [ -n "$5" ]; then
#     log_message "使用单一索引的R绘图"
#     cat "${INDEX_DIR}/sgrna.R" | sed "s/zhangt/${SAMPLE_NAME}/g" > "${WORK_DIR}/${SAMPLE_NAME}.sgrna.R"
#     chmod +x "${WORK_DIR}"/*.R
#     Rscript "${WORK_DIR}/${SAMPLE_NAME}.sgrna.R" 2>"${LOG_DIR}/${SAMPLE_NAME}.sgrna.R.log"
#     rm "${WORK_DIR}"/*.R
# else
#     log_message "使用基因和启动子的R绘图"
#     cat "${INDEX_DIR}/sgrna.gene.R" | sed "s/zhangt/${SAMPLE_NAME}/g" > "${WORK_DIR}/${SAMPLE_NAME}.sgrna.gene.R"
#     cat "${INDEX_DIR}/sgrna.promoter.R" | sed "s/zhangt/${SAMPLE_NAME}/g" > "${WORK_DIR}/${SAMPLE_NAME}.sgrna.promoter.R"
#     chmod +x "${WORK_DIR}"/*.R
#     Rscript "${WORK_DIR}/${SAMPLE_NAME}.sgrna.gene.R" 2>"${LOG_DIR}/${SAMPLE_NAME}.sgrna.gene.R.log"
#     Rscript "${WORK_DIR}/${SAMPLE_NAME}.sgrna.promoter.R" 2>"${LOG_DIR}/${SAMPLE_NAME}.sgrna.promoter.R.log"
#     rm "${WORK_DIR}"/*.R
# fi
# log_section_end "R绘图"

# 复制结果文件到结果目录
log_section_start "复制结果文件"

log_message "将结果文件复制到结果目录: ${RESULT_DIR}"
cp "${WORK_DIR}/${SAMPLE_NAME}"*.txt "${WORK_DIR}/${SAMPLE_NAME}"*.stats.txt "${RESULT_DIR}/" 2>/dev/null || {
    log_message "警告: 复制结果文件时发生错误"
}

log_section_end "复制结果文件"

# 运行完成总结
log_message "处理完成时间: $(date)"
echo "~~~~~~~处理完成~~~~~~~"

# 打印结果位置摘要
echo "结果文件位置: ${RESULT_DIR}"
echo "日志文件位置: ${LOG_DIR}"
echo "中间文件位置: ${WORK_DIR}"

# 显示处理的文件数量统计
if [ -n "$5" ]; then
    if [ -f "${WORK_DIR}/${SAMPLE_NAME}.sgRNAnumber.stats.txt" ]; then
        echo "基因统计信息:"
        cat "${WORK_DIR}/${SAMPLE_NAME}.sgRNAnumber.stats.txt"
    fi
else
    if [ -f "${WORK_DIR}/${SAMPLE_NAME}.gene.sgRNAnumber.stats.txt" ]; then
        echo "基因统计信息:"
        cat "${WORK_DIR}/${SAMPLE_NAME}.gene.sgRNAnumber.stats.txt"
    fi

    if [ -f "${WORK_DIR}/${SAMPLE_NAME}.promoter.sgRNAnumber.stats.txt" ]; then
        echo "启动子统计信息:"
        cat "${WORK_DIR}/${SAMPLE_NAME}.promoter.sgRNAnumber.stats.txt"
    fi
fi

log_message "脚本执行完成"