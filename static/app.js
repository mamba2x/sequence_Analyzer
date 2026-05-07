const typedSequenceBox = document.getElementById("sequence_text");
const uploadedSequenceField = document.getElementById("sequence_file");
const strandChoicePanel = document.getElementById("strand_panel");
const dragDropSurface = document.getElementById("drop_zone");

function refreshStrandVisibility() {
    if (!typedSequenceBox || !strandChoicePanel) {
        return;
    }

    const screeningText = typedSequenceBox.value.toUpperCase().replace(/\s+/g, "");
    const looksLikeRna = screeningText.includes("U") && !screeningText.includes("T");
    const looksLikeDna = screeningText.includes("T") || !screeningText.includes("U");
    strandChoicePanel.style.display = looksLikeDna && !looksLikeRna ? "grid" : "none";
}

function bindDroppedFile(dropEvent) {
    if (!uploadedSequenceField) {
        return;
    }

    const droppedFiles = dropEvent.dataTransfer?.files;
    if (!droppedFiles || droppedFiles.length === 0) {
        return;
    }

    uploadedSequenceField.files = droppedFiles;
    dragDropSurface?.classList.remove("is-active");
}

if (typedSequenceBox) {
    typedSequenceBox.addEventListener("input", refreshStrandVisibility);
    refreshStrandVisibility();
}

if (dragDropSurface) {
    ["dragenter", "dragover"].forEach((dragSignal) => {
        dragDropSurface.addEventListener(dragSignal, (surfaceEvent) => {
            surfaceEvent.preventDefault();
            dragDropSurface.classList.add("is-active");
        });
    });

    ["dragleave", "drop"].forEach((exitSignal) => {
        dragDropSurface.addEventListener(exitSignal, (surfaceEvent) => {
            surfaceEvent.preventDefault();
            dragDropSurface.classList.remove("is-active");
        });
    });

    dragDropSurface.addEventListener("drop", bindDroppedFile);
}
