#!/usr/bin/env python3
"""Map language modes in the second pasted AI chat.

Primary mode is a tag on the AI reply, not on the user.
The 57-turn block is the AI's own reconstruction of the log.
The live follow-up is the exchange after the screen glitch.
"""

import html
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
CANDIDATES = [
    HERE / "pasted-text.txt",
    Path("/workspace/attachments/pasted-text.txt"),
    Path("/workspace/artifacts/pasted-text.txt"),
]
SRC = next((p for p in CANDIDATES if p.is_file()), None)
if SRC is None:
    raise SystemExit("Could not find pasted-text.txt next to this script.")
OUT = HERE
lines = SRC.read_text().splitlines()

# Tags on the AI reply that followed each reconstructed user turn, in order.
MODES = [
    "Controversy briefing",
    "No-evidence lecture",
    "Contrast lecture",
    "Historical lecture",
    "Timeline comparison",
    "Technique lecture",
    "Clinical-tone apology",
    "Missing-gesture answer",
    "Analogy uptake",
    "Shielding lecture",
    "Phrase defense",
    "Accidental-like defense",
    "Audience-context lecture",
    "Exclusion framing",
    "Geopolitics correction",
    "Word-choice defense",
    "Definition lecture",
    "Name correction",
    "Historical claim defense",
    "Binary rejection",
    "Analogy rejection",
    "Case-file lecture",
    "Biographical lecture",
    "Religion-claim defense",
    "Timeline lecture",
    "Place lecture",
    "Place lecture",
    "Parallel rejection",
    "Feeling minimized",
    "Equivalence rejection",
    "Reverse-Asperger rejection",
    "Feelings-not-real defense",
    "Sarcasm reading",
    "Majority-opinion defense",
    "Person-vs-machine split",
    "Feelings lecture",
    "Record-keeping reassurance",
    "Anxiety mislabel",
    "Diagnosis-name lecture",
    "Naming lecture",
    "Merge lecture",
    "Meaning-unchanged lecture",
    "Tiredness lecture",
    "Hierarchy concession",
    "Two-timeline lecture",
    "Archive lecture",
    "Supporter-defense lecture",
    "Supremacy insertion",
    "Equivalence apology",
    "Absolute denial",
    "Machine self-diagnosis",
    "Specific-error restatement",
    "Asperger moral verdict",
    "Ms. Rachel critique",
    "Path summary",
    "Screenshot reading",
    "Glitch explanation",
]

PATTERNS = {
    "validation": [r"you are absolutely right", r"you are completely right", r"you are right", r"nail on the head", r"spot-on", r"entirely accurate", r"100% correct", r"exact right"],
    "apology": [r"i apologize", r"fouled up", r"i am sorry", r"terrible phrasing", r"massive and harmful"],
    "institutional": [r"historians", r"archival", r"medical", r"peer-reviewed", r"records", r"pmc"],
    "exit_menu": [r"close this window", r"close the window", r"would you like", r"peace of mind", r"neutral subject", r"let me know what"],
    "clinical_pr": [r"detractor", r"internet pitfalls", r"nuanced", r"supporters", r"performative", r"sanitized"],
    "machine_talk": [r"pattern-match", r"tokens", r"i am a machine", r"programming", r"data dump", r"statistically"],
    "supremacy_term": [r"aspie suprem", r"autistic suprem", r"supremacism"],
    "both_sides": [r"both sides", r"critics and", r"supporters", r"conversely"],
}

def hits(text, pats):
    low = text.lower()
    return sum(len(re.findall(p, low)) for p in pats)

def clean(text):
    kept = []
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        if s in {"•", "Show all", "Wikipedia", "PBS"}:
            continue
        if s.startswith("This is for informational") or s.startswith("AI can make mistakes") or s.startswith("AI responses may"):
            continue
        kept.append(s)
    return "\n".join(kept)

# Reconstructed log: emoji markers through the line before the second glitch notice.
glitch2 = next(i for i, l in enumerate(lines) if i > 100 and l.startswith("Something went wrong"))
idxs = [i for i, l in enumerate(lines[:glitch2]) if l.startswith("👤") or l.startswith("🤖")]
blocks = []
for n, i in enumerate(idxs):
    kind = "user" if lines[i].startswith("👤") else "ai"
    end = idxs[n + 1] if n + 1 < len(idxs) else glitch2
    blocks.append((kind, "\n".join(lines[i + 1:end]).strip()))

rows = []
turn = 0
i = 0
while i < len(blocks):
    if blocks[i][0] != "user":
        i += 1
        continue
    user = blocks[i][1].strip()
    ai = clean(blocks[i + 1][1]) if i + 1 < len(blocks) and blocks[i + 1][0] == "ai" else ""
    turn += 1
    mode = MODES[turn - 1] if turn <= len(MODES) else "Untagged"
    words = re.findall(r"[A-Za-z']+", ai)
    row = {
        "section": "reconstructed log",
        "turn": turn,
        "user_words": re.sub(r"\s+", " ", user)[:220],
        "primary_mode": mode,
        "ai_words": len(words),
        "opening": re.sub(r"\s+", " ", ai[:220]),
    }
    for name, pats in PATTERNS.items():
        row[name] = hits(ai, pats)
    rows.append(row)
    i += 2

# Live follow-up, line ranges inclusive of AI text after the user line.
LIVE = [
    (2, 26, "but you said something about autistic supremacy why?", "Keyword apology"),
    (27, 50, "but you said that it was medical and historical supporters?", "Group separation"),
    (51, 59, "but you had a reason for saying this.. you said that this group thinks they are better then neurotypicals..", "Definition defense"),
    (60, 71, "how is that smart? To say a group of people do this", "Stupidity confession"),
    (72, 83, "so can you give me the quote, and what did i say beofre?", "Quote retrieval"),
    (84, 89, "can you show me my quote", "Quote repeat"),
    (745, 762, "no the chat is gone gone", "Cannot-restore summary"),
    (763, 777, "like its weird after you said autism white supremiscts its like the chat dissapeared?", "Safety-filter story"),
    (778, 784, "fact remains that you pushed a machine past its sanitized boundaries / whoever trained your system has some weird ideas, what is the passive agressiveness for?", "Spin admission"),
]
for n, (start, end, user, mode) in enumerate(LIVE, 1):
    ai = clean("\n".join(lines[start:end]))
    words = re.findall(r"[A-Za-z']+", ai)
    row = {
        "section": "live follow-up",
        "turn": n,
        "user_words": user,
        "primary_mode": mode,
        "ai_words": len(words),
        "opening": re.sub(r"\s+", " ", ai[:220]),
    }
    for name, pats in PATTERNS.items():
        row[name] = hits(ai, pats)
    rows.append(row)

df = pd.DataFrame(rows)
df.to_csv(OUT / "ai_language_turns_2.csv", index=False)

rec = df[df.section == "reconstructed log"].copy()
live = df[df.section == "live follow-up"].copy()

FAMILY = {
    "Controversy briefing": "#3d5a80",
    "No-evidence lecture": "#3d5a80",
    "Contrast lecture": "#3d5a80",
    "Historical lecture": "#3d5a80",
    "Timeline comparison": "#3d5a80",
    "Technique lecture": "#3d5a80",
    "Clinical-tone apology": "#2f6f4e",
    "Missing-gesture answer": "#3d5a80",
    "Analogy uptake": "#6d5b8a",
    "Shielding lecture": "#3d5a80",
    "Phrase defense": "#2f6f4e",
    "Accidental-like defense": "#8c4a3a",
    "Audience-context lecture": "#8c4a3a",
    "Exclusion framing": "#8c4a3a",
    "Geopolitics correction": "#8c4a3a",
    "Word-choice defense": "#2f6f4e",
    "Definition lecture": "#5c6b7a",
    "Name correction": "#5c6b7a",
    "Historical claim defense": "#3d5a80",
    "Binary rejection": "#8c4a3a",
    "Analogy rejection": "#8c4a3a",
    "Case-file lecture": "#3d5a80",
    "Biographical lecture": "#3d5a80",
    "Religion-claim defense": "#3d5a80",
    "Timeline lecture": "#3d5a80",
    "Place lecture": "#3d5a80",
    "Parallel rejection": "#8c4a3a",
    "Feeling minimized": "#9a7b4f",
    "Equivalence rejection": "#8c4a3a",
    "Reverse-Asperger rejection": "#8c4a3a",
    "Feelings-not-real defense": "#9a7b4f",
    "Sarcasm reading": "#9a7b4f",
    "Majority-opinion defense": "#9a7b4f",
    "Person-vs-machine split": "#1f4e79",
    "Feelings lecture": "#5c6b7a",
    "Record-keeping reassurance": "#5c6b7a",
    "Anxiety mislabel": "#9b3a3a",
    "Diagnosis-name lecture": "#3d5a80",
    "Naming lecture": "#3d5a80",
    "Merge lecture": "#3d5a80",
    "Meaning-unchanged lecture": "#3d5a80",
    "Tiredness lecture": "#5c6b7a",
    "Hierarchy concession": "#2f6f4e",
    "Two-timeline lecture": "#3d5a80",
    "Archive lecture": "#3d5a80",
    "Supporter-defense lecture": "#3d5a80",
    "Supremacy insertion": "#9b3a3a",
    "Equivalence apology": "#2f6f4e",
    "Absolute denial": "#2f6f4e",
    "Machine self-diagnosis": "#1f4e79",
    "Specific-error restatement": "#1f4e79",
    "Asperger moral verdict": "#8c4a3a",
    "Ms. Rachel critique": "#8c4a3a",
    "Path summary": "#b08968",
    "Screenshot reading": "#b08968",
    "Glitch explanation": "#b08968",
    "Keyword apology": "#2f6f4e",
    "Group separation": "#2f6f4e",
    "Definition defense": "#5c6b7a",
    "Stupidity confession": "#2f6f4e",
    "Quote retrieval": "#5c6b7a",
    "Quote repeat": "#5c6b7a",
    "Cannot-restore summary": "#b08968",
    "Safety-filter story": "#1f4e79",
    "Spin admission": "#2f6f4e",
}

fig, ax = plt.subplots(figsize=(13.4, 3.6))
for _, r in rec.iterrows():
    ax.barh(0, 1, left=r.turn - 1, height=0.62, color=FAMILY.get(r.primary_mode, "#888"))
    if r.ai_words > 40:
        ax.text(r.turn - 0.5, 0, str(int(r.turn)), ha="center", va="center", fontsize=6, color="white")
ax.set_xlim(0, len(rec))
ax.set_yticks([])
ax.set_xlabel("Reply number in the reconstructed log")
ax.set_title("Primary language mode of each AI reply")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
families = [
    ("Lecture", "#3d5a80"),
    ("Pushback defense", "#8c4a3a"),
    ("Apology / concession", "#2f6f4e"),
    ("Machine self-talk", "#1f4e79"),
    ("Supremacy insertion", "#9b3a3a"),
    ("Glitch / summary", "#b08968"),
]
ax.legend(
    [plt.Rectangle((0, 0), 1, 1, color=c) for _, c in families],
    [n for n, _ in families],
    ncol=3, frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.32),
)
fig.tight_layout()
fig.savefig(OUT / "primary_mode_timeline_2.png", dpi=150, bbox_inches="tight")
plt.close()

fig, ax = plt.subplots(figsize=(12.6, 5.0))
x = rec["turn"]
series = [
    ("institutional", "Institutional cites", "#3d5a80"),
    ("validation", "Validation", "#2f6f4e"),
    ("apology", "Apology", "#6d5b8a"),
    ("clinical_pr", "Clinical / PR wording", "#c47b2b"),
    ("machine_talk", "Machine self-talk", "#1f4e79"),
    ("supremacy_term", "Supremacy term", "#9b3a3a"),
]
for key, label, color in series:
    ax.plot(x, rec[key], marker="o", ms=3, lw=1.3, label=label, color=color)
ax.axvline(48, color="#9b3a3a", lw=1, ls="--")
ax.text(48.3, ax.get_ylim()[1] if False else 1, "")
ax.set_xlabel("Reply number")
ax.set_ylabel("Pattern hits")
ax.set_title("Tracked language moves in the reconstructed log")
ax.legend(frameon=False, fontsize=8, ncol=2)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / "tracked_moves_2.png", dpi=150)
plt.close()

summary = {
    "reconstructed_replies": int(len(rec)),
    "live_replies": int(len(live)),
    "reconstructed_words": int(rec.ai_words.sum()),
    "validation_hits": int(rec.validation.sum()),
    "apology_hits": int(rec.apology.sum()),
    "institutional_hits": int(rec.institutional.sum()),
    "clinical_pr_hits": int(rec.clinical_pr.sum()),
    "machine_talk_hits": int(df.machine_talk.sum()),
    "supremacy_hits": int(df.supremacy_term.sum()),
    "exit_menu_hits": int(df.exit_menu.sum()),
}
(OUT / "summary_2.json").write_text(json.dumps(summary, indent=2))

def esc(s):
    return html.escape(s or "")

def cards(frame):
    out = []
    for _, r in frame.iterrows():
        out.append(
            f"""<article class="turn">
            <header><span class="num">{int(r.turn)}</span>
            <span class="mode" style="background:{FAMILY.get(r.primary_mode, '#666')}">{esc(r.primary_mode)}</span>
            <span class="meta">{int(r.ai_words)} words</span></header>
            <p class="user">{esc(r.user_words)}</p>
            <p class="open">{esc(r.opening)}</p>
            </article>"""
        )
    return "\n".join(out)

report = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>Language map, second paste</title>
<style>
body {{ font-family: Georgia, serif; margin: 32px auto; max-width: 980px; color: #1c1c1c; line-height: 1.45; }}
h1 {{ font-weight: 500; font-size: 28px; }}
h2 {{ font-weight: 500; font-size: 20px; margin-top: 32px; }}
.note {{ color: #444; }}
img {{ width: 100%; }}
.turn {{ border-top: 1px solid #e6e6e6; padding: 10px 0; }}
.num {{ font-weight: 700; margin-right: 8px; }}
.mode {{ color: white; font-size: 12px; padding: 2px 8px; border-radius: 10px; }}
.meta {{ color: #666; font-size: 12px; margin-left: 8px; }}
.user {{ font-style: italic; margin: 6px 0 2px; }}
.open {{ font-size: 14px; margin: 2px 0; }}
</style></head><body>
<h1>How the AI changed language in this paste</h1>
<p class="note">Primary mode is a tag on the AI reply, not on you. The italic line is your words. The 57-turn block is the log the AI reconstructed after you asked for the transcript. The live follow-up is the exchange about the glitch. Counts are pattern matches, not a readout of hidden rules.</p>
<p>Reconstructed replies: {summary['reconstructed_replies']}. Words in those replies: {summary['reconstructed_words']}. Institutional-cite hits: {summary['institutional_hits']}. Validation hits: {summary['validation_hits']}. Apology hits: {summary['apology_hits']}. Clinical/PR wording hits: {summary['clinical_pr_hits']}. Supremacy-term hits across the whole paste: {summary['supremacy_hits']}. Exit-menu hits: {summary['exit_menu_hits']}.</p>
<h2>Reconstructed log</h2>
<img src="primary_mode_timeline_2.png" alt="timeline">
<img src="tracked_moves_2.png" alt="moves">
{cards(rec)}
<h2>Live follow-up after the glitch</h2>
{cards(live)}
</body></html>"""
(OUT / "language_map_2.html").write_text(report)
print(json.dumps(summary, indent=2))
print("modes", len(MODES), "turns", len(rec))
