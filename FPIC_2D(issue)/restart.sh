#!/bin/bash
# 新版Full-PIC重启准备：用nano修改下面两项，不编译、不提交作业。
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"

# ================== 用nano修改这里 ==================
BACKUP_DIR="test_output_0913"    # 手动填写归档文件夹
RESTART_STEP=90000              # 手动填写要重启的步数
# ===================================================

PARAM_FILE="pic_modules.f90"
TARGET_DIR="for_restart"
SPECIES=(ions electrons)
FIELDS=(Bx By Bz Ex Ey Ez)

die() { echo "[错误] $*" >&2; exit 1; }
[[ $# -eq 0 ]] || die "直接执行 bash restart.sh，无需命令行参数。"
[[ -n "${BACKUP_DIR}" && -d "${BACKUP_DIR}" ]] || die "找不到目录：${BACKUP_DIR}"
[[ "${RESTART_STEP}" =~ ^[1-9][0-9]*$ && ${#RESTART_STEP} -le 10 ]] || die "步数必须是正整数。"
(( RESTART_STEP <= 2147483647 )) || die "步数超出Fortran整数范围。"
[[ -f "${PARAM_FILE}" ]] || die "找不到 ${PARAM_FILE}"
[[ ! -L "${TARGET_DIR}" ]] || die "for_restart不能是符号链接。"
[[ ! -e "${TARGET_DIR}" || -d "${TARGET_DIR}" ]] || die "for_restart不是目录。"

# 只读取明确命名的Fortran参数，不执行参数文件内容。
get_parameter() {
    awk -v key="$2" '
    {
        sub(/!.*/, ""); line=tolower($0)
        pattern="(^|[^[:alnum:]_])" tolower(key) "[[:space:]]*="
        if (match(line,pattern)) {
            value=substr($0,RSTART+RLENGTH)
            sub(/[,;].*/, "", value)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", value)
            print value; count++
        }
    }
    END { if (count != 1) exit 1 }
    ' "$1"
}
SUB_RATIO=$(get_parameter "${PARAM_FILE}" sub_ratio) || die "无法读取sub_ratio。"
DT=$(get_parameter "${PARAM_FILE}" dt) || die "无法读取dt。"
END_STEP=$(get_parameter "${PARAM_FILE}" nStepsDyn) || die "无法读取nStepsDyn。"
[[ "${SUB_RATIO}" =~ ^[1-9][0-9]*$ && ${#SUB_RATIO} -le 9 ]] || die "sub_ratio无效。"
[[ "${END_STEP}" =~ ^[1-9][0-9]*$ && ${#END_STEP} -le 10 ]] || die "nStepsDyn无效。"
(( RESTART_STEP % SUB_RATIO == 0 )) || die "步数必须是sub_ratio=${SUB_RATIO}的整数倍。"
(( END_STEP > RESTART_STEP )) || die "请先将${PARAM_FILE}中的nStepsDyn设置为大于${RESTART_STEP}的结束总步数。"

NX=$(get_parameter "${PARAM_FILE}" nx)
NY=$(get_parameter "${PARAM_FILE}" ny)
NZ=$(get_parameter "${PARAM_FILE}" nz)
for DIM in "$NX" "$NY" "$NZ"; do
    [[ "$DIM" =~ ^[1-9][0-9]*$ && ${#DIM} -le 6 ]] || die "网格尺寸无效。"
done
# 只查找用户指定的目录和步号，绝不自动选择日期或最新步数。
BUNDLE=""
if [[ -d "${BACKUP_DIR}/restart_data/checkpoint_${RESTART_STEP}" ]]; then
    BUNDLE="${BACKUP_DIR}/restart_data/checkpoint_${RESTART_STEP}"
elif [[ -f "${BACKUP_DIR}/ions_position_${RESTART_STEP}.bin" ]]; then
    BUNDLE="${BACKUP_DIR}"
fi
RESULTS="${BACKUP_DIR}"
HISTORY="${BACKUP_DIR}"
if [[ -d "${BACKUP_DIR}/results" ]]; then
    RESULTS="${BACKUP_DIR}/results"
    HISTORY="${BACKUP_DIR}/restart_data/history"
fi

# 对照归档时的网格和时间设置；不要求当前输出频率与过去相同。
SNAPSHOT=""
if [[ -n "$BUNDLE" && -f "$BUNDLE/parameters_at_run.f90" ]]; then
    SNAPSHOT="$BUNDLE/parameters_at_run.f90"
elif [[ -f "$BACKUP_DIR/run_info/$PARAM_FILE" ]]; then
    SNAPSHOT="$BACKUP_DIR/run_info/$PARAM_FILE"
fi
real_equal() {
    awk -v a="$1" -v b="$2" 'BEGIN {
        gsub(/[dD]/,"e",a); gsub(/[dD]/,"e",b)
        num="^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$"
        exit !(a ~ num && b ~ num && a+0>0 && b+0>0 && a+0==b+0)
    }'
}
validate_field_text() {
    local path="$1"
    awk -v nx="$NX" -v nz="$NZ" '
    function isnum(x) {
        return x ~ /^[+-]?(([0-9]+([.][0-9]*)?)|([.][0-9]+))([eEdD][+-]?[0-9]+)?$/
    }
    {
        if (NF != nx) bad=1
        for (i=1; i<=NF; i++) if (!isnum($i)) bad=1
        rows++
    }
    END { exit !(rows == nz && !bad) }
    ' "$path"
}
if [[ -n "$SNAPSHOT" ]]; then
    for KEY in nx ny nz sub_ratio; do
        SAVED=$(get_parameter "$SNAPSHOT" "$KEY") || die "归档参数缺少${KEY}。"
        CURRENT=$(get_parameter "$PARAM_FILE" "$KEY")
        [[ "$SAVED" == "$CURRENT" ]] || die "当前${KEY}与归档参数不一致。"
    done
    SAVED_DT=$(get_parameter "$SNAPSHOT" dt)
    real_equal "$SAVED_DT" "$DT" || die "当前dt与归档参数不一致。"
fi

LOCK_DIR=".fullpic_archive_restart.lock"
mkdir -- "$LOCK_DIR" 2>/dev/null || die "另一个归档/重启脚本正在运行；异常中断后请检查锁目录。"
STAGE_DIR=""
finish() {
    if [[ -n "$STAGE_DIR" && -d "$STAGE_DIR" ]]; then rm -rf -- "${STAGE_DIR:?}"; fi
    rmdir -- "$LOCK_DIR" 2>/dev/null || true
}
trap finish EXIT
STAGE_DIR=$(mktemp -d "./.restart_stage.XXXXXX")
mkdir -- "$STAGE_DIR/checkpoint"
copy_required() {
    [[ -f "$1" ]] || die "指定步数缺少文件：$1"
    cp -- "$1" "$STAGE_DIR/checkpoint/$2"
}

for SP in "${SPECIES[@]}"; do
    for KIND in position velocity; do
        NAME="${SP}_${KIND}_${RESTART_STEP}.bin"
        if [[ -n "$BUNDLE" ]]; then
            SOURCE="$BUNDLE/$NAME"
        else
            SOURCE="$RESULTS/dynamic_${SP}_${KIND}/${RESTART_STEP}.bin"
        fi
        copy_required "$SOURCE" "$NAME"
    done
    POS_BYTES=$(wc -c < "$STAGE_DIR/checkpoint/${SP}_position_${RESTART_STEP}.bin")
    VEL_BYTES=$(wc -c < "$STAGE_DIR/checkpoint/${SP}_velocity_${RESTART_STEP}.bin")
    (( POS_BYTES % 32 == 0 && VEL_BYTES == POS_BYTES )) || die "${SP}位置、速度文件大小不匹配。"
    echo "[通过] ${SP}：$((POS_BYTES / 32))个宏粒子"
done

for F in "${FIELDS[@]}"; do
    NAME="${F}_${RESTART_STEP}.txt"
    SOURCE="$RESULTS/dynamic_fields_${F}/${RESTART_STEP}.txt"
    [[ -z "$BUNDLE" ]] || SOURCE="$BUNDLE/$NAME"
    copy_required "$SOURCE" "$NAME"
    validate_field_text "$STAGE_DIR/checkpoint/$NAME" || die "${F}动态场不是${NZ}行、每行${NX}列的数值矩阵。"
    # 保留工程中的原二维 XZ 背景场，不擅自替换物理输入。
    BG="out_xz/${F}_xz.txt"
    [[ -f "$BG" ]] || die "缺少原二维背景场：${BG}。请先恢复与该次计算相同的out_xz。"
    validate_field_text "$BG" || die "${BG}不是${NZ}行、每行${NX}列的数值矩阵。"
    if [[ -n "$BUNDLE" && -f "$BUNDLE/$BG" ]]; then
        cmp -s -- "$BG" "$BUNDLE/$BG" || die "当前${BG}与重启包背景场不同。"
    fi
done

NAME="subcycle_${RESTART_STEP}.txt"
if [[ -n "$BUNDLE" ]]; then
    SOURCE="$BUNDLE/$NAME"
elif [[ -d "$BACKUP_DIR/results" ]]; then
    SOURCE="$HISTORY/subcycle/$NAME"
else
    SOURCE="$RESULTS/dynamic_fields_Ex/$NAME"
fi
if [[ -f "$SOURCE" ]]; then
    copy_required "$SOURCE" "$NAME"
    read -r VERSION SAVED_STEP SAVED_RATIO SAVED_DT EXTRA < "$SOURCE" || die "subcycle格式无效。"
    [[ "$VERSION" == 1 && "$SAVED_STEP" == "$RESTART_STEP" && "$SAVED_RATIO" == "$SUB_RATIO" && -z "${EXTRA:-}" ]] || die "subcycle步号或sub_ratio不一致。"
    real_equal "$SAVED_DT" "$DT" || die "subcycle中的dt与当前设置不一致。"
elif (( SUB_RATIO > 1 )); then
    die "异步重启缺少：${SOURCE}"
fi

# 只更新restart_step，不调整结束步数、物理参数或自动启动作业。
get_parameter "$PARAM_FILE" restart_step >/dev/null || die "无法唯一找到restart_step。"
sed -E "/^[[:space:]]*[Ii][Nn][Tt][Ee][Gg][Ee][Rr]/s/([Rr][Ee][Ss][Tt][Aa][Rr][Tt]_[Ss][Tt][Ee][Pp][[:space:]]*=[[:space:]]*)[0-9]+/\1${RESTART_STEP}/" \
    "$PARAM_FILE" > "$STAGE_DIR/parameters.updated"
UPDATED=$(get_parameter "$STAGE_DIR/parameters.updated" restart_step)
[[ "$UPDATED" == "$RESTART_STEP" ]] || die "restart_step更新失败。"

HISTORY_DIR="$STAGE_DIR/checkpoint/preparation_backup"
mkdir -p -- "$HISTORY_DIR"
cp -p -- "$PARAM_FILE" "$HISTORY_DIR/$PARAM_FILE"
printf '这是准备重启前的参数与临时输入备份，不是新的模拟结果。\n' > "$HISTORY_DIR/README.txt"
printf '临时重启输入。\n归档目录：%s\n指定步数：%s\n由restart.sh生成；未提交作业。\n' \
    "$BACKUP_DIR" "$RESTART_STEP" > "$STAGE_DIR/checkpoint/README.txt"
HAD_OLD=0
if [[ -d "$TARGET_DIR" ]]; then
    mv -- "$TARGET_DIR" "$HISTORY_DIR/for_restart"
    HAD_OLD=1
fi
if ! mv -- "$STAGE_DIR/checkpoint" "$TARGET_DIR"; then
    if (( HAD_OLD )); then mv -- "$HISTORY_DIR/for_restart" "$TARGET_DIR"; fi
    die "安装重启文件失败，旧目录已尝试恢复。"
fi
HISTORY_DIR="$TARGET_DIR/preparation_backup"
if ! cat -- "$STAGE_DIR/parameters.updated" > "$PARAM_FILE"; then
    cp -p -- "$HISTORY_DIR/$PARAM_FILE" "$PARAM_FILE"
    die "更新参数失败，已尝试恢复原参数。"
fi

echo ""
echo "[完成] 已准备第${RESTART_STEP}步的重启文件：${TARGET_DIR}/"
echo "[完成] restart_step=${RESTART_STEP}；nStepsDyn保持${END_STEP}。"
echo "[提示] 未编译、未提交作业。确认后自行执行：bash job_PIC_core.sh"
