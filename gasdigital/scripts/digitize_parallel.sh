#!/usr/bin/env bash
# -----------------------------------------------------------------------------
#  并行数字化：把事件号按进程数等分，每个进程只做自己那段，最后 hadd 合并。
#
#    ./gasdigital/scripts/digitize_parallel.sh pixel.root digi.root 8
#    ./gasdigital/scripts/digitize_parallel.sh pixel_height.root digi_height.root 8 --chip-y -20
#
#  每个进程用 digitize --part i/N 直接跳到自己的事件区间（靠 hits 里附带的
#  events 树，不扫描前面的数据），种子也各自错开，避免重复的随机数流。
#  合并后的 digi 树与单进程结果结构完全一致（event 号本来就是全局唯一的）。
#
#  环境变量：
#    BIN=./build/digitize   可执行文件
#    WORK=.parallel         中间目录
#    SEED=12345             第 i 个进程用 SEED+i
#    KEEP=1                 保留中间文件
# -----------------------------------------------------------------------------
set -uo pipefail

PROG=$(basename "$0")
# 仓库根目录（本脚本在 gasdigital/scripts/ 下）
ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

if [ $# -lt 3 ]; then
  cat <<EOF
用法: $PROG <hits文件> <输出digi> <进程数> [digitize 的其余参数...]

例:   $PROG pixel.root digi.root 8
      $PROG pixel_height.root digi_height.root 8 --chip-y -20 --per-cluster

环境变量:
  BIN=./build/digitize   可执行文件
  WORK=.parallel         中间目录
  SEED=12345             第 i 个进程用 SEED+i
  KEEP=1                 保留中间文件
EOF
  exit 2
fi

IN=$1
OUT=$2
JOBS=$3
shift 3
REST=("$@")
BIN=${BIN:-$ROOT_DIR/build/digitize}
WORK=${WORK:-$ROOT_DIR/.parallel}
SEED=${SEED:-12345}

if ! [[ $JOBS =~ ^[0-9]+$ ]] || [ "$JOBS" -le 0 ]; then
  echo "进程数必须是正整数：$JOBS" >&2
  exit 2
fi
if [ ! -f "$IN" ]; then
  echo "找不到输入文件：$IN" >&2
  exit 2
fi
if [ ! -x "$BIN" ]; then
  echo "找不到 $BIN（先 cmake --build build）" >&2
  exit 2
fi

mkdir -p "$WORK" || exit 1
rm -f "$WORK"/digi_part_*.root "$WORK"/digi_job_*.log

echo "[并行数字化] $IN -> $OUT，$JOBS 个进程（每进程一段事件号）"

pids=()
for ((i = 0; i < JOBS; i++)); do
  part=$(printf '%s/digi_part_%03d.root' "$WORK" "$i")
  # --seed 放在最后，覆盖用户参数里的同名选项
  "$BIN" "$IN" -o "$part" --part "$i/$JOBS" --seed "$((SEED + i))" \
    "${REST[@]}" > "$WORK/digi_job_$i.log" 2>&1 &
  pids+=($!)
done

fail=0
for idx in "${!pids[@]}"; do
  if ! wait "${pids[$idx]}"; then
    echo "[失败] worker $idx，日志尾部：" >&2
    tail -n 20 "$WORK/digi_job_$idx.log" >&2
    fail=1
  fi
done
if [ "$fail" -ne 0 ]; then
  echo "有进程失败，不合并；中间文件留在 $WORK/" >&2
  exit 1
fi

# 每个进程处理了多少事件/多少行，汇总一下便于核对
for ((i = 0; i < JOBS; i++)); do
  grep -h -E '^\[分片\]|^\[输出\]' "$WORK/digi_job_$i.log" | sed "s/^/  worker $i: /"
done

if ! command -v hadd > /dev/null 2>&1; then
  echo "找不到 hadd，中间文件在 $WORK/digi_part_*.root" >&2
  exit 1
fi

files=()
for ((i = 0; i < JOBS; i++)); do
  files+=("$(printf '%s/digi_part_%03d.root' "$WORK" "$i")")
done

echo "[合并] hadd -> $OUT"
if ! hadd -f -j "$JOBS" "$OUT" "${files[@]}" > "$WORK/digi_hadd.log" 2>&1; then
  echo "hadd 失败，见 $WORK/digi_hadd.log" >&2
  exit 1
fi

if command -v root > /dev/null 2>&1; then
  root -l -b -q -e "
    TFile f(\"$OUT\");
    auto t = (TTree*)f.Get(\"digi\");
    if (!t) { printf(\"[核对] 没有 digi 树\\n\"); return; }
    printf(\"[核对] %s: rows=%lld, 事件数=%d\\n\",
           \"$OUT\", t->GetEntries(), (int)t->GetMaximum(\"event\") + 1);
  " 2>/dev/null | grep '\[核对\]' || true
fi

if [ "${KEEP:-0}" = "1" ]; then
  echo "[完成] 中间文件保留在 $WORK/"
else
  rm -f "${files[@]}" "$WORK"/digi_job_*.log
fi
echo "[完成] $OUT"
