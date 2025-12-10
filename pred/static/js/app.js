/**
 * Southern Ground Hornbill Genetics Simulator
 * Frontend Application Logic
 */

// STATE MANAGEMENT

const state = {
    currentModel: 1,
    Ne: 500,
    generations: 50,
    results: null,
    isLoading: false
};

// Model colors for consistency
const MODEL_COLORS = {
    1: '#2563eb',  // Blue
    2: '#dc2626',  // Red
    3: '#16a34a',  // Green
    4: '#9333ea',  // Purple
    5: '#ea580c'   // Orange
};

// INITIALIZATION

document.addEventListener('DOMContentLoaded', () => {
    initializeControls();
    initializeModelSelector();
    runSimulation();
});

// CONTROL INITIALIZATION

function initializeControls() {
    // Ne slider
    const neSlider = document.getElementById('ne-slider');
    const neValue = document.getElementById('ne-value');
    
    neSlider.addEventListener('input', (e) => {
        state.Ne = parseInt(e.target.value);
        neValue.textContent = state.Ne;
    });
    
    // Generations slider
    const genSlider = document.getElementById('gen-slider');
    const genValue = document.getElementById('gen-value');
    const genUnit = document.getElementById('gen-unit');
    
    genSlider.addEventListener('input', (e) => {
        state.generations = parseInt(e.target.value);
        const years = state.generations * 26;
        genValue.textContent = state.generations;
        genUnit.textContent = `(${years} years)`;
    });
    
    // Run button
    document.getElementById('run-btn').addEventListener('click', runSimulation);
}

function initializeModelSelector() {
    const tabs = document.querySelectorAll('.model-tab');
    
    tabs.forEach(tab => {
        tab.addEventListener('click', () => {
            const modelNum = parseInt(tab.dataset.model);
            selectModel(modelNum);
        });
    });
}

function selectModel(modelNum) {
    state.currentModel = modelNum;
    
    // Update tab styling
    document.querySelectorAll('.model-tab').forEach(tab => {
        if (parseInt(tab.dataset.model) === modelNum) {
            tab.classList.add('active');
        } else {
            tab.classList.remove('active');
        }
    });
    
    // Update description
    updateModelDescription(modelNum);
    
    // Re-run simulation with new model
    runSimulation();
}

function updateModelDescription(modelNum) {
    const descriptions = {
        1: 'Baseline scenario with all wild populations (Eastern Cape, Kruger NP, KwaZulu-Natal, Limpopo)',
        2: 'Simulates loss of Eastern Cape and KwaZulu-Natal populations',
        3: 'Baseline plus 4 captive birds from PAAZA added each generation',
        4: 'Baseline plus 10 captive birds from PAAZA added each generation',
        5: 'Baseline plus 4 birds from mixed captive sources (PAAZA/AZA/EAZA) each generation'
    };
    
    document.getElementById('model-description').textContent = descriptions[modelNum];
}

// API COMMUNICATION

async function runSimulation() {
    if (state.isLoading) return;
    
    setLoading(true);
    hideError();
    
    try {
        const response = await fetch('/api/simulate', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({
                Ne: state.Ne,
                generations: state.generations,
                model: state.currentModel
            })
        });
        
        if (!response.ok) {
            const error = await response.json();
            throw new Error(error.error || 'Simulation failed');
        }
        
        state.results = await response.json();
        displayResults(state.results);
        
    } catch (error) {
        showError(error.message);
    } finally {
        setLoading(false);
    }
}

// UI STATE MANAGEMENT

function setLoading(isLoading) {
    state.isLoading = isLoading;
    
    const loadingEl = document.getElementById('loading');
    const resultsEl = document.getElementById('results');
    const runBtn = document.getElementById('run-btn');
    
    if (isLoading) {
        loadingEl.classList.remove('hidden');
        resultsEl.classList.add('hidden');
        runBtn.disabled = true;
    } else {
        loadingEl.classList.add('hidden');
        runBtn.disabled = false;
    }
}

function showError(message) {
    const errorEl = document.getElementById('error');
    errorEl.textContent = message;
    errorEl.classList.remove('hidden');
}

function hideError() {
    document.getElementById('error').classList.add('hidden');
}

// RESULTS DISPLAY

function displayResults(results) {
    document.getElementById('results').classList.remove('hidden');
    
    displayStatistics(results);
    plotCharts(results);
}

function displayStatistics(results) {
    const initial = results.initial;
    const data = results.Ho; // Array of Ho values
    const finalHo = data[data.length - 1];
    const finalHe = results.He[results.He.length - 1];
    const finalF = results.F[results.F.length - 1];
    const finalNa = results.Na[results.Na.length - 1];
    const finalN = results.population[results.population.length - 1];
    
    // Calculate changes
    const hoChange = ((finalHo - initial.Ho) / initial.Ho * 100).toFixed(1);
    const naChange = ((finalNa - initial.Na) / initial.Na * 100).toFixed(1);
    const nChange = ((finalN - initial.N) / initial.N * 100).toFixed(1);
    
    const statsHtml = `
        <div class="stat-card">
            <div class="stat-label">Heterozygosity (Ho)</div>
            <div class="stat-value">${finalHo.toFixed(4)}</div>
            <div class="stat-change ${hoChange < 0 ? 'negative' : 'positive'}">
                ${hoChange}% from initial (${initial.Ho.toFixed(4)})
            </div>
        </div>
        
        <div class="stat-card">
            <div class="stat-label">Inbreeding Coefficient (F)</div>
            <div class="stat-value">${finalF.toFixed(4)}</div>
            <div class="stat-change ${finalF > 0.25 ? 'warning' : ''}">
                ${finalF > 0.25 ? '⚠️ Above concern threshold (0.25)' : '✓ Below concern threshold'}
            </div>
        </div>
        
        <div class="stat-card">
            <div class="stat-label">Allelic Richness (Na)</div>
            <div class="stat-value">${finalNa.toFixed(2)}</div>
            <div class="stat-change ${naChange < 0 ? 'negative' : 'positive'}">
                ${naChange}% from initial (${initial.Na.toFixed(2)})
            </div>
        </div>
        
        <div class="stat-card">
            <div class="stat-label">Population Size</div>
            <div class="stat-value">${Math.round(finalN)}</div>
            <div class="stat-change">
                ${nChange}% from initial (${initial.N})
            </div>
        </div>
    `;
    
    document.getElementById('stats-grid').innerHTML = statsHtml;
}

// CHART PLOTTING

function plotCharts(results) {
    plotHeterozygosity(results);
    plotInbreeding(results);
    plotAlleles(results);
    plotPopulation(results);
}

function plotHeterozygosity(results) {
    const trace = {
        x: results.generations,
        y: results.Ho,
        mode: 'lines',
        name: 'Observed (Ho)',
        line: {
            color: MODEL_COLORS[state.currentModel],
            width: 3
        }
    };
    
    const layout = {
        title: {
            text: 'Heterozygosity Over Time',
            font: { size: 18, color: '#1f2937' }
        },
        xaxis: {
            title: 'Generation',
            gridcolor: '#e5e7eb'
        },
        yaxis: {
            title: 'Heterozygosity',
            range: [0, 0.6],
            gridcolor: '#e5e7eb'
        },
        plot_bgcolor: '#ffffff',
        paper_bgcolor: '#ffffff',
        hovermode: 'closest'
    };
    
    Plotly.newPlot('chart-heterozygosity', [trace], layout, {responsive: true});
}

function plotInbreeding(results) {
    const trace = {
        x: results.generations,
        y: results.F,
        mode: 'lines',
        name: 'Inbreeding (F)',
        line: {
            color: MODEL_COLORS[state.currentModel],
            width: 3
        }
    };
    
    const layout = {
        title: {
            text: 'Inbreeding Coefficient Over Time',
            font: { size: 18, color: '#1f2937' }
        },
        xaxis: {
            title: 'Generation',
            gridcolor: '#e5e7eb'
        },
        yaxis: {
            title: 'Inbreeding Coefficient (F)',
            range: [0, 1],
            gridcolor: '#e5e7eb'
        },
        plot_bgcolor: '#ffffff',
        paper_bgcolor: '#ffffff',
        shapes: [{
            type: 'line',
            x0: 0,
            x1: results.generations[results.generations.length - 1],
            y0: 0.25,
            y1: 0.25,
            line: {
                color: '#f59e0b',
                width: 2,
                dash: 'dash'
            }
        }],
        annotations: [{
            x: results.generations[Math.floor(results.generations.length / 2)],
            y: 0.25,
            text: 'Concern Threshold (F = 0.25)',
            showarrow: false,
            yshift: 10,
            font: { color: '#f59e0b' }
        }]
    };
    
    Plotly.newPlot('chart-inbreeding', [trace], layout, {responsive: true});
}

function plotAlleles(results) {
    const trace = {
        x: results.generations,
        y: results.Na,
        mode: 'lines',
        name: 'Alleles (Na)',
        line: {
            color: MODEL_COLORS[state.currentModel],
            width: 3
        }
    };
    
    const layout = {
        title: {
            text: 'Allelic Richness Over Time',
            font: { size: 18, color: '#1f2937' }
        },
        xaxis: {
            title: 'Generation',
            gridcolor: '#e5e7eb'
        },
        yaxis: {
            title: 'Mean Number of Alleles',
            range: [0, 8],
            gridcolor: '#e5e7eb'
        },
        plot_bgcolor: '#ffffff',
        paper_bgcolor: '#ffffff'
    };
    
    Plotly.newPlot('chart-alleles', [trace], layout, {responsive: true});
}

function plotPopulation(results) {
    const trace = {
        x: results.generations,
        y: results.population,
        mode: 'lines',
        name: 'Population',
        line: {
            color: MODEL_COLORS[state.currentModel],
            width: 3
        }
    };
    
    const layout = {
        title: {
            text: 'Population Size Projection',
            font: { size: 18, color: '#1f2937' }
        },
        xaxis: {
            title: 'Generation',
            gridcolor: '#e5e7eb'
        },
        yaxis: {
            title: 'Population Size (N)',
            gridcolor: '#e5e7eb'
        },
        plot_bgcolor: '#ffffff',
        paper_bgcolor: '#ffffff'
    };
    
    Plotly.newPlot('chart-population', [trace], layout, {responsive: true});
}