document.addEventListener("DOMContentLoaded", function () {
    const searchInput = document.getElementById("searchInput");
    const geneFilter = document.getElementById("geneFilter");
    const tableBody = document.getElementById("trials-body");

    let allStudies = [];
    let geneChartInstance = null;

    // Fetch data from Flask API
    fetch("/api/studies")
        .then(response => {
            if (!response.ok) throw new Error("Network response error");
            return response.json();
        })
        .then(data => {
            allStudies = data;
            renderTable(allStudies);
            renderChart(allStudies);
        })
        .catch(error => {
            console.error("Fetch error:", error);
            tableBody.innerHTML = `<tr><td colspan="8" class="text-danger text-center">Failed to load data from server.</td></tr>`;
        });

    function renderTable(studies) {
        tableBody.innerHTML = "";

        if (studies.length === 0) {
            tableBody.innerHTML = `<tr><td colspan="8" class="text-center">No matching records found.</td></tr>`;
            return;
        }

        studies.forEach(study => {
            const geneBadge = study.primary_gene !== 'None' 
                ? `<span class="badge bg-primary">${study.primary_gene}</span>`
                : `<span class="text-muted">None</span>`;

            const patentBadge = study.patent_count !== 'N/A'
                ? `<span class="badge bg-info text-dark">${study.patent_count} patents</span>`
                : `<span class="text-muted">N/A</span>`;

            let indexBadge = `<span class="text-muted">N/A</span>`;
            if (study.ip_index_label.includes('Saturated')) {
                indexBadge = `<span class="badge bg-warning text-dark">${study.ip_index_label}</span>`;
            } else if (study.ip_index_label.includes('Unmet Need')) {
                indexBadge = `<span class="badge bg-danger">${study.ip_index_label}</span>`;
            }

            const row = document.createElement("tr");
            row.className = "trial-row";
            row.innerHTML = `
                <td><strong>${study.nct_id}</strong></td>
                <td class="text-start">${study.title}</td>
                <td><span class="badge bg-success">${study.status}</span></td>
                <td>${geneBadge}</td>
                <td>${patentBadge}</td>
                <td>${indexBadge}</td>
                <td><small>${study.conditions}</small></td>
                <td><small>${study.sponsor}</small></td>
            `;
            tableBody.appendChild(row);
        });
    }

    function renderChart(studies) {
        const geneMetrics = {};

        // Aggregate patent and trial counts per target gene
        studies.forEach(s => {
            const gene = s.primary_gene;
            if (gene && gene !== "None") {
                if (!geneMetrics[gene]) {
                    geneMetrics[gene] = {
                        trials: 0,
                        patents: typeof s.patent_count === "number" ? s.patent_count : 0
                    };
                }
                geneMetrics[gene].trials += 1;
            }
        });

        const labels = Object.keys(geneMetrics);
        const trialCounts = labels.map(g => geneMetrics[g].trials);
        const patentCounts = labels.map(g => geneMetrics[g].patents);

        const ctx = document.getElementById("geneChart").getContext("2d");

        // Destroy old instance if re-rendering
        if (geneChartInstance) {
            geneChartInstance.destroy();
        }

        geneChartInstance = new Chart(ctx, {
            type: "bar",
            data: {
                labels: labels,
                datasets: [
                    {
                        label: "Total Patents (Europe PMC)",
                        data: patentCounts,
                        backgroundColor: "rgba(13, 202, 240, 0.7)",
                        borderColor: "rgba(13, 202, 240, 1)",
                        borderWidth: 1,
                        yAxisID: "yPatents"
                    },
                    {
                        label: "Active Clinical Trials",
                        data: trialCounts,
                        backgroundColor: "rgba(13, 110, 253, 0.7)",
                        borderColor: "rgba(13, 110, 253, 1)",
                        borderWidth: 1,
                        yAxisID: "yTrials"
                    }
                ]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    yPatents: {
                        type: "linear",
                        position: "left",
                        title: { display: true, text: "Patent Count" },
                        beginAtZero: true
                    },
                    yTrials: {
                        type: "linear",
                        position: "right",
                        title: { display: true, text: "Trial Count" },
                        beginAtZero: true,
                        grid: { drawOnChartArea: false }
                    }
                }
            }
        });
    }

    function filterTable() {
        const query = searchInput.value.toLowerCase();
        const selectedGene = geneFilter.value.toLowerCase();

        const filtered = allStudies.filter(study => {
            const fullText = `${study.nct_id} ${study.title} ${study.primary_gene} ${study.sponsor}`.toLowerCase();
            const matchesSearch = fullText.includes(query);
            const matchesGene = selectedGene === "" || study.primary_gene.toLowerCase() === selectedGene;

            return matchesSearch && matchesGene;
        });

        renderTable(filtered);
        renderChart(filtered);
    }

    searchInput.addEventListener("keyup", filterTable);
    geneFilter.addEventListener("change", filterTable);
});