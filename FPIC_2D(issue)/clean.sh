#!/bin/bash
# 作业结束后执行：按日期归档全部动态输出，不运行模拟
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
shopt -s nullglob

die() { echo "[错误] $*" >&2; exit 1; }
[[ $# -eq 0 ]] || die "用法：bash clean.sh"
OUTPUTS=(dynamic_* diagnostics)
PRESENT=()
for ITEM in "${OUTPUTS[@]}"; do
    [[ -e "$ITEM" ]] || continue
    [[ -d "$ITEM" && ! -L "$ITEM" ]] || die "输出路径不是普通目录：${ITEM}"
    PRESENT+=("$ITEM")
done
(( ${#PRESENT[@]} > 0 )) || die "没有输出目录可归档，未改动任何文件。"
for ITEM in for_restart restart_history; do
    [[ ! -L "$ITEM" ]] || die "临时目录不能是符号链接：${ITEM}"
    [[ ! -e "$ITEM" || -d "$ITEM" ]] || die "临时路径不是目录：${ITEM}"
done
LOCK_DIR=".fullpic_archive_restart.lock"
mkdir -- "$LOCK_DIR" 2>/dev/null || die "归档或重启脚本正在运行；异常中断后请检查锁目录。"
STAGE=""
COMMITTED=0
FROM=()
TO=()
finish() {
    CODE=$?
    trap - EXIT
    if (( ! COMMITTED )); then
        for (( I=${#FROM[@]}-1; I>=0; I-- )); do
            if [[ -e "${TO[I]}" ]]; then
                mv -- "${TO[I]}" "${FROM[I]}" || {
                    echo "[错误] 恢复失败，请保留并检查暂存目录：${STAGE}" >&2
                    rmdir -- "$LOCK_DIR" 2>/dev/null || true
                    exit 1
                }
            fi
        done
        [[ -z "$STAGE" || ! -d "$STAGE" ]] || rm -rf -- "${STAGE:?}"
    fi
    rmdir -- "$LOCK_DIR" 2>/dev/null || true
    exit "$CODE"
}
trap finish EXIT
move_item() {
    mkdir -p -- "$(dirname -- "$2")"
    mv -- "$1" "$2"
    FROM+=("$1"); TO+=("$2")
}
ARCHIVE="test_output_$(date +%m%d)"
if [[ -e "$ARCHIVE" ]]; then ARCHIVE="${ARCHIVE}_$(date +%H%M%S)_$$"; fi
[[ ! -e "$ARCHIVE" ]] || die "归档目录已存在：${ARCHIVE}"
STAGE=$(mktemp -d "./.archive_stage.XXXXXX")
mkdir -p -- "$STAGE/results" "$STAGE/restart_data/history" "$STAGE/restart_temp" "$STAGE/run_info"

# 保存本轮参数与源代码，便于以后核对网格、dt和sub_ratio。
for ITEM in *.f90 job_PIC_core.sh clean.sh cleanup.sh restart.sh; do
    [[ ! -f "$ITEM" ]] || cp -p -- "$ITEM" "$STAGE/run_info/"
done
for ITEM in "${PRESENT[@]}"; do
    move_item "$ITEM" "$STAGE/results/$ITEM"
done
# subcycle是异步续算所需的元数据。
for ITEM in "$STAGE/results/dynamic_fields_Ex"/subcycle_*.txt; do
    move_item "$ITEM" "$STAGE/restart_data/history/subcycle/$(basename -- "$ITEM")"
done
# 所有可删除的旧重启输入和准备备份，集中到同一个文件夹。
for ITEM in for_restart restart_history; do
    [[ ! -e "$ITEM" ]] || move_item "$ITEM" "$STAGE/restart_temp/$ITEM"
done
for ITEM in PIC_CORE_Zhang.log PIC_CORE_Zhang.err PIC_CORE_Zhang_qsub_job.sh .latest_fullpic_archive; do
    [[ ! -f "$ITEM" ]] || move_item "$ITEM" "$STAGE/run_info/$ITEM"
done
cat > "$STAGE/README.txt" <<'EOF'
这是本轮模拟数据归档。
results/：位置、速度、六个动态场分量和诊断数据，用于分析。
restart_data/：各输出步的subcycle等原始续算元数据。需要以后续算时请保留。
restart_temp/：已经使用过的for_restart及准备备份，不是本轮新结果，可整文件夹删除。
run_info/：本轮源码、参数和日志。

续算：nano restart.sh，填写本归档目录和所需步数，再执行bash restart.sh。
所选步必须同时具有完整的粒子位置、速度和六个场输出；异步还需subcycle。
restart.sh会生成新的临时for_restart，不会编译或提交作业。
原二维背景场out_xz仍保留在工程目录，续算时必须与本轮背景场一致。
本归档不再重复复制一套最新检查点；可以选择任意具备完整输出的历史步。
EOF
printf '这个文件夹仅包含已经使用过的重启临时输入和准备备份，可以整体删除。\n不要误删相邻的restart_data，它是以后续算必需的数据。\n' > "$STAGE/restart_temp/README.txt"
printf '这里保存subcycle历史。它们不是可丢弃的临时文件。\n重启时还需results中同一步的粒子位置、速度和六个动态场，以及原二维背景场out_xz。\n' > "$STAGE/restart_data/README.txt"
mv -- "$STAGE" "$ARCHIVE"
COMMITTED=1

# 编译产物和测试生成器的临时可执行文件无需随模拟数据保存。
for ITEM in *.mod *.smod *.o *.obj \
    PIC_CORE_Zhang PIC_CORE_Zhang.exe \
    field_harris_xz field_harris_xz.exe \
    particle_harris_xz particle_harris_xz.exe; do
    [[ ! -f "$ITEM" ]] || rm -f -- "$ITEM"
done
echo "[完成] 本轮输出已归档到 ${ARCHIVE}/"
echo "[提示] ${ARCHIVE}/restart_temp/ 可以整体删除。"
echo "[下一步] nano restart.sh → bash restart.sh → bash job_PIC_core.sh"
