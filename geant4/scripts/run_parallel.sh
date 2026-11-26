#!/usr/bin/env bash
# -----------------------------------------------------------------------------
#  并行跑 pixel：把 N 个事件拆给若干独立进程，各写各的文件，最后 hadd 合并。
#
#    ./geant4/scripts/run_parallel.sh geant4/macros/run_height.mac 20000 15 pixel_height.root
#
#  每个进程用不同的事件号偏移（第 i 个从前面所有进程的事件数之和起），
#  所以 hadd 之后 event 全局唯一、整体递增，plot.py / plot_digi.py /
#  digitize 的"按事件号分组"全部照旧可用。
#
#  分支名、分支类型、每行的物理量都不变：脚本只改宏里的输出文件名、
#  事件号偏移、随机数种子和 /run/beamOn 的事件数。
#
#  环境变量：
#    BIN=./build/pixel    可执行文件路径
#    WORK=.parallel       中间文件目录（每个进程一个宏 + 一个 root + 一个日志）
#    KEEP=1               合并后保留中间文件（默认删掉）
# -----------------------------------------------------------------------------
set -uo pipefail

PROG=$(basename "$0")
# 仓库根目录（本脚本在 geant4/scripts/ 下），默认路径都相对它，从哪运行都行
ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)

if [ $# -lt 4 ]; then
  cat <<EOF
用法: $PROG <基础宏> <总事件数> <进程数> <输出文件>

例:   $PROG geant4/macros/run.mac 5000 8 pixel.root
      $PROG geant4/macros/run_height.mac 20000 15 pixel_height.root

环境变量:
  BIN=./build/pixel   可执行文件
  WORK=.parallel      中间文件目录
  KEEP=1              保留中间文件
EOF
  exit 2
fi

BASE=$1
TOTAL=$2
JOBS=$3
OUT=$4
BIN=${BIN:-$ROOT_DIR/build/pixel}
WORK=${WORK:-$ROOT_DIR/.parallel}

# ---- 参数检查 ----
for v in "$TOTAL" "$JOBS"; do
  if ! [[ $v =~ ^[0-9]+$ ]]; then
    echo "事件数和进程数必须是正整数：$v" >&2
    exit 2
  fi
done
if [ "$TOTAL" -le 0 ] || [ "$JOBS" -le 0 ]; then
  echo "事件数和进程数必须 > 0" >&2
  exit 2
fi
if [ ! -f "$BASE" ]; then
  echo "找不到基础宏：$BASE" >&2
  exit 2
fi
if [ ! -x "$BIN" ]; then
  echo "找不到可执行文件：$BIN（先 cmake --build build）" >&2
  exit 2
fi
if [ "$JOBS" -gt "$TOTAL" ]; then
  echo "[警告] 进程数 $JOBS > 事件数 $TOTAL，按 $TOTAL 个进程跑" >&2
  JOBS=$TOTAL
fi

mkdir -p "$WORK" || exit 1
rm -f "$WORK"/part_*.root "$WORK"/job_*.mac "$WORK"/job_*.log

# ---- 基础宏里的随机种子（没有就用默认值）----
read -r S1 S2 <<< "$(awk '/^[[:space:]]*\/random\/setSeeds/ {print $2, $3; exit}' "$BASE")"
S1=${S1:-12345}
S2=${S2:-67890}

# ---- 基础宏去掉"每个进程要自己定"的几行，其余原样保留 ----
#      留下的顺序不变，所以 /run/initialize 仍在这些命令前面。
FILTERED="$WORK/base.mac"
grep -Ev '^[[:space:]]*(/output/fileName|/output/eventOffset|/random/setSeeds|/run/beamOn)' \
  "$BASE" > "$FILTERED" || true

# ---- 分片：尽量平均，余数从前面的进程开始各多 1 个 ----
BASE_N=$(( TOTAL / JOBS ))
REM=$(( TOTAL % JOBS ))

echo "[并行] 基础宏 $BASE -> $TOTAL 个事件 / $JOBS 个进程（每进程 $BASE_N~$((BASE_N + 1)) 个）"

pids=()
offset=0
for ((i = 0; i < JOBS; i++)); do
  n=$BASE_N
  if [ "$i" -lt "$REM" ]; then n=$((n + 1)); fi

  mac="$WORK/job_$i.mac"
  part=$(printf '%s/part_%03d.root' "$WORK" "$i")

  cat "$FILTERED" > "$mac"
  {
    echo ""
    echo "# ---- 并行分片 $i/$JOBS（run_parallel.sh 生成，别手改）----"
    echo "/output/fileName  $part"
    echo "/output/eventOffset $offset"
    # 每个进程换一对种子，避免跑的是同一串随机数
    echo "/random/setSeeds $((S1 + i)) $((S2 + i))"
    echo "/run/beamOn $n"
  } >> "$mac"

  "$BIN" "$mac" > "$WORK/job_$i.log" 2>&1 &
  pids+=($!)

  echo "  worker $i: event $offset..$((offset + n - 1))  ->  $part"
  offset=$((offset + n))
done

# ---- 等所有进程结束；任何一个失败就不合并 ----
fail=0
for idx in "${!pids[@]}"; do
  if ! wait "${pids[$idx]}"; then
    echo "[失败] worker $idx，日志尾部：" >&2
    tail -n 20 "$WORK/job_$idx.log" >&2
    fail=1
  fi
done
if [ "$fail" -ne 0 ]; then
  echo "有进程失败，不合并；中间文件留在 $WORK/" >&2
  exit 1
fi

# ---- 合并：按 worker 编号顺序，事件号因此整体递增 ----
if ! command -v hadd > /dev/null 2>&1; then
  echo "找不到 hadd（ROOT 的 bin 目录不在 PATH 里？）" >&2
  echo "中间文件在 $WORK/part_*.root" >&2
  exit 1
fi

files=()
for ((i = 0; i < JOBS; i++)); do
  files+=("$(printf '%s/part_%03d.root' "$WORK" "$i")")
done

echo "[合并] hadd -> $OUT"
if ! hadd -f -j "$JOBS" "$OUT" "${files[@]}" > "$WORK/hadd.log" 2>&1; then
  echo "hadd 失败，见 $WORK/hadd.log" >&2
  exit 1
fi

# ---- 核对：条数与事件号最大值（有 root 才算，没有就跳过）----
if command -v root > /dev/null 2>&1; then
  root -l -b -q -e "
    TFile f(\"$OUT\");
    auto t = (TTree*)f.Get(\"hits\");
    if (!t) { printf(\"[核对] 没有 hits 树\\n\"); return; }
    printf(\"[核对] %s: entries=%lld, 事件数=%d\\n\",
           \"$OUT\", t->GetEntries(), (int)t->GetMaximum(\"event\") + 1);
  " 2>/dev/null | grep '\[核对\]' || true
fi

if [ "${KEEP:-0}" = "1" ]; then
  echo "[完成] 中间文件保留在 $WORK/"
else
  rm -f "${files[@]}" "$WORK"/job_*.mac "$WORK"/job_*.log "$FILTERED"
  echo "[完成] 中间文件已清理（KEEP=1 可保留）"
fi
echo "[完成] $OUT"
