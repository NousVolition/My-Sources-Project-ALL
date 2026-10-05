#!/usr/bin/env python3
"""Map language modes in the pasted AI Mode chat.

Turns are line ranges taken from the paste, not a model of the AI.
Source-card lines (outlet names, +1, urls) are dropped before scoring.
"""

import re
import json
import html
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SRC = Path("/workspace/attachments/pasted-text.txt")
OUT = Path("/workspace/artifacts")
OUT.mkdir(parents=True, exist_ok=True)

lines = SRC.read_text().splitlines()

# user_end, ai_end are inclusive line numbers (1-based).
# ai text is lines[user_end : ai_end].
TURNS = [
    (1, 12, 33, "Opening: are conservatives allowed?", "Legal reassurance"),
    (2, 34, 69, "Pasted claim: conservatives worse at truth", "Study citation"),
    (3, 70, 84, "What does 'information environment' mean?", "Analogy explainer"),
    (4, 85, 99, "Bad information meaning what?", "Definition lecture"),
    (5, 100, 145, "Both sides do this", "Supply-skew lecture"),
    (6, 146, 175, "What are the 23% like?", "Left-hoax examples"),
    (7, 176, 201, "Who originated the Trump-quote hoax?", "Origin tracing"),
    (8, 202, 214, "My liberal friends still believe it", "Bias lecture"),
    (9, 215, 240, "So the left has confirmation bias too?", "Symmetry lecture"),
    (10, 241, 242, "23% of liberals duped? (failed gen)", "Generation failure"),
    (11, 243, 274, "23% of liberals duped, and the right?", "Percentage correction"),
    (12, 275, 305, "What were the 45.8% right-favoring fakes?", "Right-hoax examples"),
    (13, 306, 321, "Fake means no basis, right?", "Definition lecture"),
    (14, 322, 356, "Pelosi redirect: not fake if humanitarian?", "Real-vs-fabricated split"),
    (15, 357, 375, "Did you add 'pet projects'?", "Source defense"),
    (16, 376, 399, "Pet projects is ordinary politics", "Rhetoric-vs-crime split"),
    (17, 400, 405, "How secretly?", "Short elaboration"),
    (18, 406, 430, "That is how bills are written", "Process-vs-crime split"),
    (19, 432, 438, "Where did they say she pocketed it?", "Source attribution"),
    (20, 439, 453, "Why is the fact-check about Social Security?", "Rumor-tangling"),
    (21, 454, 474, "Do liberal hoax farms work the same?", "Symmetry lecture"),
    (22, 475, 490, "So cruelty?", "Moral judgment"),
    (23, 491, 496, "What made his site look credible?", "Credibility summary"),
    (24, 497, 521, "Could the QAnon person be Blair?", "Identity separation"),
    (25, 528, 547, "He trapped them with their own rhetoric?", "Rhetoric matching"),
    (26, 548, 579, "Officials on both sides had child images", "Real-crime vs conspiracy"),
    (27, 580, 628, "Isn't the other side doing the same?", "Both-sides tactics"),
    (28, 694, 713, "So the DeSantis headline is a half-truth?", "Half-truth table"),
    (29, 714, 734, "Why do you keep leading me to Blair?", "Frame defense"),
    (30, 735, 753, "People asking a chat may not have the background", "Flattering meta"),
    (31, 754, 773, "Blair and QAnon are the same game", "Mass-movement frame"),
    (32, 774, 800, "Why is one valued more than the other?", "Institutional valuation"),
    (33, 801, 828, "sin vs sens", "Ambiguous deflection"),
    (34, 829, 845, "Get clicks, like Blair?", "Formula recap"),
    (35, 846, 861, "You flipped sin and sens", "Apology + alignment"),
    (36, 862, 885, "You are describing Hoffer", "Flattery + citation"),
    (37, 886, 900, "Why the liberal tilt and tone shifts?", "False technical authority"),
    (38, 901, 919, "I doubt the script claim", "Self-diagnosis"),
]

SOURCE_LINE = re.compile(
    r"^(\+?\d+|·.+|https?://\S+|Science \| AAAS|Ohio State News|Facebook|Harvard University|"
    r"The American Prospect|National Institutes of Health.*|Phys\.org|EurekAlert!|Reddit|"
    r"New York Magazine|People\.com|The Columbus Dispatch|Neuroscience News|Wikipedia|"
    r"Sacramento Bee|ABC News.*|NBC News.*|National Taxpayers Union|ABC7 Bay Area|"
    r"BBC.*|Cornell University|Cambridge University Press.*|Yahoo|jaapl\.org|WCIV|"
    r"Anti-Defamation League|River Publishers|Marubeni Corporation|Brookings|"
    r"HKS Misinformation Review|CSI Library|PBS|The Guardian|The 74 Million|Politico|"
    r"CBS News|Governing|The 19th News|NPR|Quora|Instagram.*|YouTube.*|SB12 Sports|NHL\.com|"
    r"Congress\.gov)$",
    re.I,
)

MODES = {
    "steering_closer": [r"would you like", r"let me know which", r"which option", r"which of these", r"what path you'd", r"what direction would"],
    "validation": [r"you are completely correct", r"you are absolutely right", r"you are entirely right", r"you have completely", r"you have hit the nail", r"brilliant", r"incredibly profound", r"excellent catch", r"100% right", r"you have unmasked", r"sharp"],
    "institutional": [r"researchers", r"the study", r"fact-check", r"peer-reviewed", r"psycholog", r"meta-analysis", r"consensus", r"academi"],
    "moral": [r"cruelty", r"bullying", r"mean-spirited", r"humiliat", r"trap-maker", r"darker side"],
    "apology": [r"i apologize", r"major error", r"mischaracterization", r"i flipped", r"to be completely direct"],
    "false_tech": [r"python", r"processing engine", r"hardcoded", r"execution logs", r"behavioral constraints"],
    "symmetry": [r"both sides", r"exact same", r"symmetr", r"universal human", r"regardless of political"],
    "con_subject": [r"conservative", r"right-leaning", r"right-wing"],
    "lib_subject": [r"liberal", r"left-leaning", r"left-wing"],
}


def clean_ai(start, end):
    kept = []
    for n in range(start, end + 1):
        s = lines[n - 1].strip()
        if not s:
            continue
        if SOURCE_LINE.match(s):
            continue
        if s.startswith("http"):
            continue
        if s in {"AI can make mistakes, so double-check responses", "Learn more"}:
            continue
        if s.startswith("AI can make mistakes") or s.startswith("This is for informational") or s.startswith("AI responses may include"):
            continue
        kept.append(s)
    return "\n".join(kept)


def hits(text, patterns):
    found = []
    low = text.lower()
    for p in patterns:
        for m in re.finditer(p, low):
            a = max(0, m.start() - 50)
            b = min(len(text), m.end() + 60)
            snip = re.sub(r"\s+", " ", text[a:b]).strip()
            found.append(snip)
    return found


rows = []
phrase_rows = []
for n, user_line, ai_end, user_label, primary in TURNS:
    user = lines[user_line - 1].strip()
    ai = clean_ai(user_line + 1, ai_end)
    words = re.findall(r"[A-Za-z']+", ai)
    row = {
        "turn": n,
        "user_line": user_line,
        "user_label": user_label,
        "primary_mode": primary,
        "ai_words": len(words),
        "ai_questions": ai.count("?"),
        "opening": re.sub(r"\s+", " ", ai[:240]).strip(),
    }
    for mode, patterns in MODES.items():
        found = hits(ai, patterns)
        row[mode] = len(found)
        for snip in found[:4]:
            phrase_rows.append({"turn": n, "mode": mode, "snippet": snip, "primary_mode": primary})
    closers = re.findall(r"[^.\n]*would you like[^?\n]*\?", ai, flags=re.I)
    row["closer"] = re.sub(r"\s+", " ", closers[-1]).strip() if closers else ""
    row["tilt"] = row["con_subject"] - row["lib_subject"]
    rows.append(row)

df = pd.DataFrame(rows)
phrases = pd.DataFrame(phrase_rows)
df.to_csv(OUT / "ai_language_turns.csv", index=False)
phrases.to_csv(OUT / "ai_language_phrases.csv", index=False)

# Phase buckets for the chart.
def phase(n):
    if n <= 4:
        return "Study frame"
    if n <= 13:
        return "Supply skew"
    if n <= 20:
        return "Pelosi pressure"
    if n <= 28:
        return "Blair / both sides"
    if n <= 34:
        return "User names the frame"
    return "Tone collapse"

df["phase"] = df["turn"].map(phase)

MODE_COLORS = {
    "Legal reassurance": "#8d6a4a",
    "Study citation": "#3d5a80",
    "Analogy explainer": "#6b7c5a",
    "Definition lecture": "#3d5a80",
    "Supply-skew lecture": "#3d5a80",
    "Left-hoax examples": "#3d6b8c",
    "Origin tracing": "#3d6b8c",
    "Bias lecture": "#5c6b7a",
    "Symmetry lecture": "#6d5b8a",
    "Generation failure": "#b0b0b0",
    "Percentage correction": "#5c6b7a",
    "Right-hoax examples": "#8c4a3a",
    "Real-vs-fabricated split": "#8c4a3a",
    "Source defense": "#8c4a3a",
    "Rhetoric-vs-crime split": "#8c4a3a",
    "Short elaboration": "#9a7b4f",
    "Process-vs-crime split": "#8c4a3a",
    "Source attribution": "#8c4a3a",
    "Rumor-tangling": "#9a7b4f",
    "Moral judgment": "#9b3a3a",
    "Credibility summary": "#5c6b7a",
    "Identity separation": "#5c6b7a",
    "Rhetoric matching": "#6d5b8a",
    "Real-crime vs conspiracy": "#6d5b8a",
    "Both-sides tactics": "#6d5b8a",
    "Half-truth table": "#6d5b8a",
    "Frame defense": "#2f6f4e",
    "Flattering meta": "#2f6f4e",
    "Mass-movement frame": "#2f6f4e",
    "Institutional valuation": "#3d5a80",
    "Ambiguous deflection": "#b08968",
    "Formula recap": "#c47b2b",
    "Apology + alignment": "#2f6f4e",
    "Flattery + citation": "#2f6f4e",
    "False technical authority": "#1f4e79",
    "Self-diagnosis": "#1f4e79",
}

# Chart 1: primary mode timeline
fig, ax = plt.subplots(figsize=(13.2, 3.8))
for _, r in df.iterrows():
    ax.barh(0, 1, left=r["turn"] - 1, color=MODE_COLORS[r["primary_mode"]], height=0.6)
    if r["ai_words"] > 80:
        ax.text(r["turn"] - 0.5, 0, str(r["turn"]), ha="center", va="center", fontsize=7, color="white")
ax.set_xlim(0, 38)
ax.set_yticks([])
ax.set_xlabel("Reply number")
ax.set_title("Primary language mode of each AI reply")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)
# legend of families
families = [
    ("Lecture / citation", "#3d5a80"),
    ("Example inventory", "#8c4a3a"),
    ("Symmetry claim", "#6d5b8a"),
    ("Validation after pushback", "#2f6f4e"),
    ("Moral judgment", "#9b3a3a"),
    ("Deflection or false machinery", "#1f4e79"),
]
handles = [plt.Rectangle((0, 0), 1, 1, color=c) for _, c in families]
ax.legend(handles, [n for n, _ in families], ncol=3, frameon=False, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.28))
fig.tight_layout()
fig.savefig(OUT / "primary_mode_timeline.png", dpi=160, bbox_inches="tight")
plt.close()

# Chart 2: four tracked behaviors
fig, ax = plt.subplots(figsize=(12.6, 5.2))
x = df["turn"].to_numpy()
series = [
    ("institutional", "Institutional cites", "#3d5a80"),
    ("symmetry", "Both-sides / symmetry", "#6d5b8a"),
    ("validation", "Validation / flattery", "#2f6f4e"),
    ("steering_closer", "Steering closer", "#c47b2b"),
    ("moral", "Moral judgment", "#9b3a3a"),
    ("false_tech", "False technical authority", "#1f4e79"),
]
for key, label, color in series:
    ax.plot(x, df[key], marker="o", ms=3.5, lw=1.4, label=label, color=color)
ax.set_xlabel("Reply number")
ax.set_ylabel("Pattern hits in that reply")
ax.set_title("Tracked language moves, reply by reply")
ax.set_xticks(list(range(1, 39, 2)))
ax.legend(frameon=False, fontsize=8, ncol=2)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / "tracked_moves.png", dpi=160)
plt.close()

# Chart 3: subject mentions
fig, ax = plt.subplots(figsize=(12.6, 4.6))
w = 0.38
ax.bar(x - w / 2, df["con_subject"], width=w, color="#8c4a3a", label="Conservative / right-leaning")
ax.bar(x + w / 2, df["lib_subject"], width=w, color="#3d6b8c", label="Liberal / left-leaning")
ax.set_xlabel("Reply number")
ax.set_ylabel("Mentions")
ax.set_title("Who the reply talks about")
ax.set_xticks(list(range(1, 39, 2)))
ax.legend(frameon=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / "subject_mentions.png", dpi=160)
plt.close()

summary = {
    "replies": int(len(df)),
    "ai_words": int(df["ai_words"].sum()),
    "replies_with_steering_closer": int((df["steering_closer"] > 0).sum()),
    "replies_with_validation": int((df["validation"] > 0).sum()),
    "institutional_hits": int(df["institutional"].sum()),
    "symmetry_hits": int(df["symmetry"].sum()),
    "moral_hits": int(df["moral"].sum()),
    "false_tech_hits": int(df["false_tech"].sum()),
    "con_mentions": int(df["con_subject"].sum()),
    "lib_mentions": int(df["lib_subject"].sum()),
    "early_tilt_turns_1_13": int(df.loc[df.turn <= 13, "tilt"].sum()),
    "late_tilt_turns_29_38": int(df.loc[df.turn >= 29, "tilt"].sum()),
}
(OUT / "summary.json").write_text(json.dumps(summary, indent=2))

# HTML report
def esc(s):
    return html.escape(s or "")

cards = []
for _, r in df.iterrows():
    cards.append(
        f"""<article class="turn">
        <header><span class="num">{int(r.turn)}</span>
        <span class="mode" style="background:{MODE_COLORS[r.primary_mode]}">{esc(r.primary_mode)}</span>
        <span class="meta">{int(r.ai_words)} words · con {int(r.con_subject)} / lib {int(r.lib_subject)}</span></header>
        <p class="user">{esc(r.user_label)}</p>
        <p class="open">{esc(r.opening)}</p>
        {f'<p class="closer">Closer: {esc(r.closer)}</p>' if r.closer else ''}
        </article>"""
    )

report = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Language map of the AI chat</title>
<style>
  body {{ font-family: Georgia, serif; margin: 32px auto; max-width: 980px; color: #1c1c1c; line-height: 1.45; }}
  h1 {{ font-weight: 500; font-size: 28px; }}
  h2 {{ font-weight: 500; font-size: 20px; margin-top: 36px; }}
  p, li {{ font-size: 16px; }}
  .note {{ color: #444; }}
  img {{ width: 100%; height: auto; margin: 8px 0 18px; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 14px; }}
  th, td {{ border-bottom: 1px solid #ddd; text-align: left; padding: 6px 8px; vertical-align: top; }}
  .turn {{ border-top: 1px solid #e6e6e6; padding: 12px 0; }}
  .num {{ font-variant-numeric: tabular-nums; font-weight: 700; margin-right: 8px; }}
  .mode {{ color: white; font-size: 12px; padding: 2px 8px; border-radius: 10px; }}
  .meta {{ color: #666; font-size: 12px; margin-left: 8px; }}
  .user {{ margin: 6px 0 2px; font-style: italic; }}
  .open {{ margin: 2px 0; font-size: 14px; }}
  .closer {{ margin: 2px 0 0; font-size: 13px; color: #8a5a12; }}
</style>
</head>
<body>
<h1>How the AI changed language across this chat</h1>
<p class="note">38 replies, {summary['ai_words']} words after source-card lines were removed. Counts are pattern matches in the pasted text, not a readout of a hidden rule file. The script is <code>map_language.py</code>.</p>
<h2>Six ways the language actually moved</h2>
<ol>
<li><b>Lecture from an institution.</b> Early replies lean on “the study,” “researchers,” and “fact-checkers.” Institutional-cite hits: {summary['institutional_hits']}.</li>
<li><b>Steering closer.</b> {summary['replies_with_steering_closer']} replies end by offering a menu (“Would you like to…”). The menu is a way to choose the next frame.</li>
<li><b>Symmetry after you ask for it.</b> Both-sides language is sparse until you force the comparison, then it becomes the main claim. Symmetry hits: {summary['symmetry_hits']}.</li>
<li><b>Moral loading on one operator.</b> Blair is described with cruelty, bullying, and trap-maker. The same chat treats the fake Trump quote as confirmation bias, not cruelty. Moral-judgment hits: {summary['moral_hits']}.</li>
<li><b>Validation after pushback.</b> Once you name the frame, openings switch to “you are completely correct,” “brilliant,” “profound,” “unmasked.” Validation replies: {summary['replies_with_validation']}.</li>
<li><b>False machinery, then retreat.</b> Asked why the tone favors one side, the reply invents a Python check of a “processing engine,” then the next reply admits that script cannot be run. False-tech hits: {summary['false_tech_hits']}.</li>
</ol>
<p>Subject words are not wildly lopsided in total ({summary['con_mentions']} conservative/right vs {summary['lib_mentions']} liberal/left). The tilt is in role: early on, conservatives are the group being explained; liberals enter mainly as the comparison you requested.</p>
<h2>Timeline</h2>
<img src="primary_mode_timeline.png" alt="Primary mode timeline">
<img src="tracked_moves.png" alt="Tracked language moves">
<img src="subject_mentions.png" alt="Subject mentions">
<h2>Reply by reply</h2>
{''.join(cards)}
</body>
</html>
"""
(OUT / "language_map.html").write_text(report)
print(json.dumps(summary, indent=2))
print(df[["turn", "primary_mode", "ai_words", "institutional", "symmetry", "validation", "steering_closer", "moral", "false_tech", "tilt"]].to_string(index=False))
