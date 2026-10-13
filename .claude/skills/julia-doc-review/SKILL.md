---
name: julia-doc-review
description: Use after a change of the README, the manual or public docstrings of this Julia package, to review the text as a reader who has no access to the source.
---

# Julia documentation review

The review is independent of the writing: the reviewer gets the text and reader tasks, not the source and not the reasons of the author. Use the `julia-doc-reader` agent for that isolation. The usual outcome of a review is to cut, clarify or move text; adding text needs a reader question that the text left open.

## Reader questions

Pose three to five questions that the intended reader has, and for a task page add one task to complete from the page alone.

- What is the package, and when is it the right tool?
- How do I install it and run the first example?
- What does this call return, in which units, and what does it throw?
- What is the main limit I must know before I rely on it?
- Where is the detail for <topic>?

Mark each answer: **supported** (quote the passage), **ambiguous** (quote, say what is unclear) or **absent**. Do not answer from knowledge of the repository.

## Per-section inspection

| Check | Action when it fails |
|---|---|
| The section answers a reader question or enables an action | delete it |
| A term is not standard Julia or domain vocabulary and is not defined | define it in one sentence, or replace it |
| A sentence states no fact | delete it |
| The same explanation exists elsewhere | keep one home, link from the other |
| Internal detail in a user page | move it to `CONTRACT.md` or delete it |
| An example does not run or shows invented output | make it run, or make it a plain block without output |
| A number in prose | move it to a table |

## Verdict

Report the questions with their marks, the list of edits made (cut, clarified, moved), and what stays open. Length, the count of headings and apparent completeness are never reasons to approve.
