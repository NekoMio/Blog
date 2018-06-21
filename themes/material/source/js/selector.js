renderWithKaTeX = function(e) {
    renderMathInElement(
        e, {
            delimiters: [
                { left: "$$", right: "$$", display: true },
                { left: "$", right: "$", display: false }
            ]
        }
    );
}

document.addEventListener("DOMContentLoaded", function() {
    renderWithKaTeX(document.body);
});
