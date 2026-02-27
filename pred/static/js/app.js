/**
 * Southern Ground Hornbill Genetics Simulator
 * Frontend Application Logic
 */

// STATE MANAGEMENT

const state = {
    currentModel: 1,
    Ne: 500,
    generations: 30,
    lambda: 1.00,
    stochastic: false,
    monteCarlo: false,
    mcReplicates: 100,
    extinctionThreshold: 100,
    results: null,
    isLoading: false
};

// Active SSE connection for Monte Carlo streaming
let mcEventSource = null;

// Model colors for consistency
const MODEL_COLORS = {
    1: '#2563eb',  // Blue
    2: '#dc2626',  // Red
    3: '#16a34a',  // Green
    4: '#9333ea',  // Purple
    5: '#ea580c',  // Orange
    6: '#0d9488'   // Teal
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

    // Lambda (growth rate) slider
    const lambdaSlider = document.getElementById('lambda-slider');
    const lambdaValue = document.getElementById('lambda-value');
    const lambdaUnit = document.getElementById('lambda-unit');

    lambdaSlider.addEventListener('input', (e) => {
        state.lambda = parseFloat(e.target.value);
        lambdaValue.textContent = state.lambda.toFixed(2);

        // Update descriptive label
        if (state.lambda < 0.98) {
            lambdaUnit.textContent = '(declining)';
        } else if (state.lambda > 1.02) {
            lambdaUnit.textContent = '(growing)';
        } else {
            lambdaUnit.textContent = '(stable)';
        }
    });

    // Stochastic toggle
    const stochasticToggle = document.getElementById('stochastic-toggle');
    const stochasticInfo = document.getElementById('stochastic-info');

    stochasticToggle.addEventListener('change', (e) => {
        state.stochastic = e.target.checked;
        stochasticInfo.style.display = e.target.checked ? 'block' : 'none';
    });

    // Monte Carlo toggle
    const mcToggle = document.getElementById('mc-toggle');
    const mcOptions = document.getElementById('mc-options');

    mcToggle.addEventListener('change', (e) => {
        state.monteCarlo = e.target.checked;
        mcOptions.style.display = e.target.checked ? 'block' : 'none';
        // MC mode implies stochastic; disable/enable the stochastic toggle accordingly
        stochasticToggle.disabled = e.target.checked;
        stochasticToggle.parentElement.style.opacity = e.target.checked ? '0.5' : '1';
    });

    // Monte Carlo replicates slider
    const mcReplicatesSlider = document.getElementById('mc-replicates-slider');
    const mcReplicatesValue = document.getElementById('mc-replicates-value');

    mcReplicatesSlider.addEventListener('input', (e) => {
        state.mcReplicates = parseInt(e.target.value);
        mcReplicatesValue.textContent = state.mcReplicates;
    });

    // Quasi-extinction threshold input
    const mcThreshold = document.getElementById('mc-threshold');

    mcThreshold.addEventListener('change', (e) => {
        const val = parseInt(e.target.value);
        if (val >= 1 && val <= 10000) {
            state.extinctionThreshold = val;
        }
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
        3: 'Baseline plus 4 South African captive birds (SA captive zoos) added each generation',
        4: 'Baseline plus 10 South African captive birds (SA captive zoos) added each generation',
        5: 'Baseline plus 4 birds from mixed captive sources (South African + USA + European zoos) each generation',
        6: 'Mixed SA sourcing: 6 captive + 2 KZN wild + 2 EC wild birds per generation'
    };

    document.getElementById('model-description').textContent = descriptions[modelNum];
}

// API COMMUNICATION

async function runSimulation() {
    if (state.isLoading) return;

    // Cancel any in-progress MC stream before starting a new one
    if (mcEventSource) {
        mcEventSource.close();
        mcEventSource = null;
    }

    setLoading(true);
    hideError();

    if (state.monteCarlo) {
        runMCStream();   // SSE path — manages its own loading state
        return;
    }

    try {
        const response = await fetch('/api/simulate', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                Ne: state.Ne,
                generations: state.generations,
                lambda: state.lambda,
                model: state.currentModel,
                stochastic: state.stochastic
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

function runMCStream() {
    const params = new URLSearchParams({
        model: state.currentModel,
        Ne: state.Ne,
        generations: state.generations,
        lambda: state.lambda,
        n_replicates: state.mcReplicates,
        extinction_threshold: state.extinctionThreshold,
    });

    setMCProgress(0, state.mcReplicates);

    mcEventSource = new EventSource(`/api/monte-carlo/stream?${params}`);

    mcEventSource.onmessage = (event) => {
        let data;
        try { data = JSON.parse(event.data); }
        catch { return; }

        if (data.type === 'error') {
            mcEventSource.close();
            mcEventSource = null;
            setLoading(false);
            hideMCProgress();
            showError(data.message || 'Monte Carlo stream error');
            return;
        }

        // Update progress bar
        setMCProgress(data.completed, data.total);

        // Render partial or final results
        state.results = data;
        displayMCResults(data);

        if (data.type === 'complete') {
            mcEventSource.close();
            mcEventSource = null;
            setLoading(false);
            hideMCProgress();
        }
    };

    mcEventSource.onerror = () => {
        // Only show error if we haven't already received a complete event
        if (mcEventSource) {
            mcEventSource.close();
            mcEventSource = null;
            setLoading(false);
            hideMCProgress();
            if (!state.results) {
                showError('Monte Carlo stream disconnected. Try reducing replicates.');
            }
        }
    };
}

// UI STATE MANAGEMENT

function setLoading(isLoading) {
    state.isLoading = isLoading;

    const loadingEl = document.getElementById('loading');
    const resultsEl = document.getElementById('results');
    const runBtn = document.getElementById('run-btn');

    if (isLoading) {
        document.getElementById('loading-text').textContent = 'Running simulation';
        loadingEl.classList.remove('hidden');
        resultsEl.classList.add('hidden');
        runBtn.disabled = true;
    } else {
        loadingEl.classList.add('hidden');
        runBtn.disabled = false;
    }
}

function setMCProgress(completed, total) {
    const wrap = document.getElementById('mc-progress-bar-wrap');
    const bar  = document.getElementById('mc-progress-bar');
    const text = document.getElementById('loading-text');
    wrap.style.display = 'block';
    const pct = total > 0 ? Math.round(completed / total * 100) : 0;
    bar.style.width = pct + '%';
    text.textContent = `Monte Carlo: ${completed} / ${total} replicates (${pct}%)`;
}

function hideMCProgress() {
    document.getElementById('mc-progress-bar-wrap').style.display = 'none';
    document.getElementById('mc-progress-bar').style.width = '0%';
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
    document.getElementById('chart-extinction-container').style.display = 'none';

    displayStatistics(results);
    plotCharts(results);
}

function displayMCResults(results) {
    document.getElementById('results').classList.remove('hidden');
    document.getElementById('chart-extinction-container').style.display = 'block';

    displayMCStatistics(results);
    // Use Plotly.react() for smooth progressive updates (avoids full re-layout flicker)
    plotMCCharts(results);
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

function displayMCStatistics(results) {
    const initial = results.initial;
    const m = results.metrics;
    const extProb = (results.extinction_probability_final * 100).toFixed(1);
    const extColor = results.extinction_probability_final > 0.1 ? '#dc2626' :
                     results.extinction_probability_final > 0.01 ? '#f59e0b' : '#16a34a';

    const finalHoMedian = m.Ho.p50[m.Ho.p50.length - 1];
    const finalNaMedian = m.Na.p50[m.Na.p50.length - 1];
    const finalNMedian  = m.population.p50[m.population.p50.length - 1];
    const finalFMedian  = m.F.p50[m.F.p50.length - 1];

    const hoChange = ((finalHoMedian - initial.Ho) / initial.Ho * 100).toFixed(1);
    const naChange = ((finalNaMedian - initial.Na) / initial.Na * 100).toFixed(1);
    const nChange  = ((finalNMedian  - initial.N)  / initial.N  * 100).toFixed(1);

    const statsHtml = `
        <div class="stat-card">
            <div class="stat-label">Heterozygosity (Ho) — median</div>
            <div class="stat-value">${finalHoMedian.toFixed(4)}</div>
            <div class="stat-change ${hoChange < 0 ? 'negative' : 'positive'}">
                ${hoChange}% from initial (${initial.Ho.toFixed(4)})
            </div>
        </div>

        <div class="stat-card">
            <div class="stat-label">Inbreeding Coefficient (F) — median</div>
            <div class="stat-value">${finalFMedian.toFixed(4)}</div>
            <div class="stat-change ${finalFMedian > 0.25 ? 'warning' : ''}">
                ${finalFMedian > 0.25 ? '⚠️ Above concern threshold (0.25)' : '✓ Below concern threshold'}
            </div>
        </div>

        <div class="stat-card">
            <div class="stat-label">Allelic Richness (Na) — median</div>
            <div class="stat-value">${finalNaMedian.toFixed(2)}</div>
            <div class="stat-change ${naChange < 0 ? 'negative' : 'positive'}">
                ${naChange}% from initial (${initial.Na.toFixed(2)})
            </div>
        </div>

        <div class="stat-card">
            <div class="stat-label">Population Size — median</div>
            <div class="stat-value">${Math.round(finalNMedian)}</div>
            <div class="stat-change">${nChange}% from initial (${initial.N})</div>
        </div>

        <div class="stat-card" style="border-left: 4px solid ${extColor};">
            <div class="stat-label">Quasi-Extinction Probability</div>
            <div class="stat-value" style="color: ${extColor};">${extProb}%</div>
            <div class="stat-change">
                N &lt; ${results.extinction_threshold} at generation ${results.generations[results.generations.length - 1]}
                &nbsp;(${results.n_replicates} replicates)
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
    
    Plotly.react('chart-heterozygosity', [trace], layout, {responsive: true});
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
    
    Plotly.react('chart-inbreeding', [trace], layout, {responsive: true});
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
    
    Plotly.react('chart-alleles', [trace], layout, {responsive: true});
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

    Plotly.react('chart-population', [trace], layout, {responsive: true});
}

// MONTE CARLO CHART HELPERS

/**
 * Build Plotly traces for a confidence band chart.
 * Uses fill: 'tonexty' to shade between percentile pairs.
 * Order: upper-p95, lower-p5 (fills outer band), upper-p75, lower-p25 (fills inner band), median, mean.
 */
function buildMCTraces(gens, metrics, color) {
    const alpha = (a) => {
        // Convert hex color to rgba
        const r = parseInt(color.slice(1,3), 16);
        const g = parseInt(color.slice(3,5), 16);
        const b = parseInt(color.slice(5,7), 16);
        return `rgba(${r},${g},${b},${a})`;
    };

    return [
        // Outer band upper boundary (p95) — rendered first
        {
            x: gens, y: metrics.p95,
            mode: 'lines', line: { width: 0, color: color },
            showlegend: false, hoverinfo: 'skip', name: 'p95'
        },
        // Outer band lower boundary (p5) — fills to p95
        {
            x: gens, y: metrics.p5,
            mode: 'lines', fill: 'tonexty', fillcolor: alpha(0.12),
            line: { width: 0, color: color },
            showlegend: false, hoverinfo: 'skip', name: '90% range'
        },
        // Inner band upper boundary (p75)
        {
            x: gens, y: metrics.p75,
            mode: 'lines', line: { width: 0, color: color },
            showlegend: false, hoverinfo: 'skip', name: 'p75'
        },
        // Inner band lower boundary (p25) — fills to p75
        {
            x: gens, y: metrics.p25,
            mode: 'lines', fill: 'tonexty', fillcolor: alpha(0.28),
            line: { width: 0, color: color },
            showlegend: false, hoverinfo: 'skip', name: '50% range'
        },
        // Median line
        {
            x: gens, y: metrics.p50,
            mode: 'lines', line: { color: color, width: 2.5 },
            name: 'Median'
        },
        // Mean line (dashed)
        {
            x: gens, y: metrics.mean,
            mode: 'lines', line: { color: color, width: 1.5, dash: 'dot' },
            name: 'Mean', opacity: 0.7
        }
    ];
}

function plotMCCharts(results) {
    const gens = results.generations;
    const color = MODEL_COLORS[state.currentModel];
    const subtitle = `(${results.n_replicates} replicates, shading = 50%/90% CI)`;

    // Heterozygosity
    Plotly.react('chart-heterozygosity',
        buildMCTraces(gens, results.metrics.Ho, color),
        {
            title: { text: `Heterozygosity Over Time<br><sub>${subtitle}</sub>`, font: { size: 18, color: '#1f2937' } },
            xaxis: { title: 'Generation', gridcolor: '#e5e7eb' },
            yaxis: { title: 'Heterozygosity', range: [0, 0.7], gridcolor: '#e5e7eb' },
            plot_bgcolor: '#ffffff', paper_bgcolor: '#ffffff', hovermode: 'x unified'
        },
        { responsive: true }
    );

    // Inbreeding
    const lastGen = gens[gens.length - 1];
    Plotly.react('chart-inbreeding',
        buildMCTraces(gens, results.metrics.F, color),
        {
            title: { text: `Inbreeding Coefficient Over Time<br><sub>${subtitle}</sub>`, font: { size: 18, color: '#1f2937' } },
            xaxis: { title: 'Generation', gridcolor: '#e5e7eb' },
            yaxis: { title: 'Inbreeding Coefficient (F)', range: [0, 1], gridcolor: '#e5e7eb' },
            plot_bgcolor: '#ffffff', paper_bgcolor: '#ffffff', hovermode: 'x unified',
            shapes: [{ type: 'line', x0: 0, x1: lastGen, y0: 0.25, y1: 0.25,
                        line: { color: '#f59e0b', width: 2, dash: 'dash' } }],
            annotations: [{ x: gens[Math.floor(gens.length / 2)], y: 0.25,
                            text: 'Concern threshold (F = 0.25)', showarrow: false,
                            yshift: 10, font: { color: '#f59e0b' } }]
        },
        { responsive: true }
    );

    // Allelic richness
    Plotly.react('chart-alleles',
        buildMCTraces(gens, results.metrics.Na, color),
        {
            title: { text: `Allelic Richness Over Time<br><sub>${subtitle}</sub>`, font: { size: 18, color: '#1f2937' } },
            xaxis: { title: 'Generation', gridcolor: '#e5e7eb' },
            yaxis: { title: 'Mean Number of Alleles', range: [0, 10], gridcolor: '#e5e7eb' },
            plot_bgcolor: '#ffffff', paper_bgcolor: '#ffffff', hovermode: 'x unified'
        },
        { responsive: true }
    );

    // Population size
    Plotly.react('chart-population',
        buildMCTraces(gens, results.metrics.population, color),
        {
            title: { text: `Population Size Projection<br><sub>${subtitle}</sub>`, font: { size: 18, color: '#1f2937' } },
            xaxis: { title: 'Generation', gridcolor: '#e5e7eb' },
            yaxis: { title: 'Population Size (N)', gridcolor: '#e5e7eb' },
            plot_bgcolor: '#ffffff', paper_bgcolor: '#ffffff', hovermode: 'x unified'
        },
        { responsive: true }
    );

    // Extinction probability
    const extProbPct = results.extinction_probability.map(p => p * 100);
    Plotly.react('chart-extinction',
        [
            {
                x: gens, y: extProbPct,
                mode: 'lines', fill: 'tozeroy',
                fillcolor: 'rgba(220,38,38,0.12)',
                line: { color: '#dc2626', width: 2.5 },
                name: `P(N < ${results.extinction_threshold})`
            }
        ],
        {
            title: {
                text: `Quasi-Extinction Probability (N < ${results.extinction_threshold} birds)<br><sub>${results.n_replicates} replicates</sub>`,
                font: { size: 18, color: '#1f2937' }
            },
            xaxis: { title: 'Generation', gridcolor: '#e5e7eb' },
            yaxis: { title: 'Probability (%)', range: [0, 100], gridcolor: '#e5e7eb' },
            plot_bgcolor: '#ffffff', paper_bgcolor: '#ffffff',
            shapes: [{
                type: 'line', x0: 0, x1: lastGen, y0: 10, y1: 10,
                line: { color: '#f59e0b', width: 2, dash: 'dash' }
            }],
            annotations: [{
                x: gens[Math.floor(gens.length / 2)], y: 10,
                text: '10% alert threshold', showarrow: false,
                yshift: 10, font: { color: '#f59e0b' }
            }]
        },
        { responsive: true }
    );
}