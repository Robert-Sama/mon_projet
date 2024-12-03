document.getElementById("file-form").addEventListener("submit", async (e) => {
    e.preventDefault();

    const formData = new FormData();
    formData.append("query_file", document.getElementById("query_file").files[0]);
    formData.append("db_file", document.getElementById("db_file").files[0]);
    formData.append("g", document.getElementById("g").value);
    formData.append("E", document.getElementById("E").value);
    formData.append("ss", document.getElementById("ss").value);

    const resultsDiv = document.getElementById("results");
    resultsDiv.innerHTML = "Processing...";

    try {
        const response = await fetch("/process", {
            method: "POST",
            body: formData,
        });

        if (response.ok) {
            const data = await response.json();

            // Afficher le tableau HTML retourné par le serveur
            resultsDiv.innerHTML = data.table;
        } else {
            const error = await response.json();
            resultsDiv.innerHTML = `<p style="color: red;">Erreur : ${error.error}</p>`;
        }
    } catch (err) {
        console.error("Erreur lors du traitement :", err);
        resultsDiv.innerHTML = `<p style="color: red;">Erreur lors du traitement des données.</p>`;
    }
});
