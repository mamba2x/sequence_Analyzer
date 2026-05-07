# Project 2 Assignment

## Run the app

```powershell
python app.py
```

Then open `http://127.0.0.1:5000` in your browser.

## Deploy on Render

- Push this folder to GitHub.
- Create a new `Web Service` on Render and connect the repository.
- Render can use the included `render.yaml`, or you can set:
  - Build command: `pip install -r requirements.txt`
  - Start command: `gunicorn app:biology_portal`

## What it does

- Accepts DNA or RNA input by typing, file upload, and drag-and-drop.
- Detects DNA, RNA, or invalid input with a plain-English explanation.
- Supports DNA coding strand and template strand transcription.
- Translates mRNA into codons and amino acids.
- Displays the amino-acid chain with full names and abbreviations.
- Characterizes the resulting protein and attempts a live UniProt lookup.

## Notes

- FASTA headers beginning with `>` are ignored automatically.
- If UniProt cannot be reached, the app still works and explains that the live lookup was unavailable.
- Variable names in the Python source were written to be distinct and descriptive instead of repetitive.
