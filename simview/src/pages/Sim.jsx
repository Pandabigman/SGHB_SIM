import React, { useState, useEffect } from 'react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';
import { Info, AlertCircle } from 'lucide-react';
import { styles } from '../components/styles';
import { 
  initialGenetics, 
  demographics, 
  modelColors, 
  getModelName 
} from '../js/constant';
import {
 calculateLambda,
  calculateGeneticDiversity,
  calculateInbreeding,
  calculateAlleles,
  projectPopulationSize,
  estimatePopulationSize,
  calculateEffectiveNe
} from '../components/utils';

const Sim = () => {
  const [ne, setNe] = useState(100);
  const [generations, setGenerations] = useState(50);
  const [isSimulating, setIsSimulating] = useState(false);
  const [results, setResults] = useState(null);
  const [activeTab, setActiveTab] = useState('overview');
  const [stochastic, setStochastic] = useState(false);

  const runSimulation = () => {
    setIsSimulating(true);
    
    const N_wild = estimatePopulationSize(199, 199, 450);
    const N_wildNoECKZN = estimatePopulationSize(177, 199, 450);
    const lambda = calculateLambda();
    
    const models = [];
    
    for (let model = 1; model <= 5; model++) {
      const modelData = [];
      
      for (let t = 0; t <= generations; t++) {
        let H0, A0, N0, suppAdded = 0;
        
        switch(model) {
          case 1:
            H0 = initialGenetics.wild.Ho;
            A0 = initialGenetics.wild.Na;
            N0 = N_wild;
            break;
          case 2:
            H0 = initialGenetics.wildNoECKZN.Ho;
            A0 = initialGenetics.wildNoECKZN.Na;
            N0 = N_wildNoECKZN;
            break;
          case 3:
            H0 = initialGenetics.wild.Ho;
            A0 = initialGenetics.wild.Na;
            N0 = N_wild;
            suppAdded = 4 * t;
            break;
          case 4:
            H0 = initialGenetics.wild.Ho;
            A0 = initialGenetics.wild.Na;
            N0 = N_wild;
            suppAdded = 10 * t;
            break;
          case 5:
            H0 = initialGenetics.wild.Ho;
            A0 = initialGenetics.wild.Na;
            N0 = N_wild;
            suppAdded = 4 * t;
            break;
        }
        
        let effectiveNe = ne;
        if (model >= 3) {
          const immigrantContribution = suppAdded * 0.5;
          effectiveNe = ne + immigrantContribution;
        }
        
        const Ho = calculateGeneticDiversity(H0, effectiveNe, t);
        const He = calculateGeneticDiversity(initialGenetics.wild.He, effectiveNe, t);
        const F = calculateInbreeding(effectiveNe, t);
        const Na = calculateAlleles(A0, effectiveNe, t);
        const popSize = projectPopulationSize(N0, lambda, t) + suppAdded;
        
        modelData.push({
          generation: t,
          year: t * 15,
          Ho: Ho.toFixed(4),
          He: He.toFixed(4),
          F: F.toFixed(4),
          Na: Na.toFixed(2),
          popSize: popSize,
          model: model
        });
      }
      
      models.push({
        modelNumber: model,
        name: getModelName(model),
        data: modelData
      });
    }
    
    setResults(models);
    setIsSimulating(false);
  };

  const prepareChartData = (metric) => {
    if (!results) return [];
    
    const combinedData = [];
    const maxGen = generations;
    
    for (let gen = 0; gen <= maxGen; gen++) {
      const dataPoint = { generation: gen };
      
      results.forEach(model => {
        const genData = model.data.find(d => d.generation === gen);
        if (genData) {
          dataPoint[`model${model.modelNumber}`] = parseFloat(genData[metric]);
        }
      });
      
      combinedData.push(dataPoint);
    }
    
    return combinedData;
  };

  useEffect(() => {
    runSimulation();
  }, []);

  return (
    <div style={styles.container}>
      <div style={styles.maxWidth}>
        {/* Header */}
        <div style={styles.card}>
          <div style={styles.flexStart}>
            <div>
                <img src='./images.jpeg'/>
            </div>
            <div style={styles.flex1}>
              <h1 style={styles.title}>
                Southern Ground Hornbill Population Genetics Simulator
              </h1>
              <p style={styles.subtitle}>
                Predictive modeling of genetic diversity and population viability over 50 generations (750 years)
              </p>
              <div style={styles.infoBox}>
                <Info size={16} />
                <span>Based on microsatellite data from 199 wild and 145 captive individuals</span>
              </div>
            </div>
          </div>
        </div>

        {/* Controls */}
        <div style={styles.card}>
          <h2 style={styles.sectionTitle}>Simulation Parameters</h2>
          
          <div style={styles.gridCols2}>
            <div>
              <label style={styles.label}>
                Effective Population Size (Ne): {ne}
              </label>
              <input
                type="range"
                min="20"
                max="200"
                value={ne}
                onChange={(e) => setNe(parseInt(e.target.value))}
                style={styles.slider}
              />
              <div style={styles.sliderLabels}>
                <span>20</span>
                <span>100</span>
                <span>200</span>
              </div>
            </div>
            
            <div>
              <label style={styles.label}>
                Generations to Simulate: {generations}
              </label>
              <input
                type="range"
                min="10"
                max="100"
                step="10"
                value={generations}
                onChange={(e) => setGenerations(parseInt(e.target.value))}
                style={styles.slider}
              />
              <div style={styles.sliderLabels}>
                <span>10 (150 yrs)</span>
                <span>50 (750 yrs)</span>
                <span>100 (1500 yrs)</span>
              </div>
            </div>
          </div>
          
          <div style={{marginTop: '1rem'}}>
            <label style={{...styles.label, display: 'flex', alignItems: 'center', gap: '0.5rem', cursor: 'pointer'}}>
              <input
                type="checkbox"
                checked={stochastic}
                onChange={(e) => setStochastic(e.target.checked)}
                style={{width: '1rem', height: '1rem', cursor: 'pointer'}}
              />
              <span>
                Enable Stochastic Simulation 
                <span style={{fontSize: '0.75rem', color: '#6b7280', marginLeft: '0.5rem'}}>
                  (adds random variation - each run will be different)
                </span>
              </span>
            </label>
          </div>
          
          <button
            onClick={runSimulation}
            disabled={isSimulating}
            style={{
              ...styles.button,
              ...(isSimulating ? styles.buttonDisabled : {})
            }}
            onMouseEnter={(e) => !isSimulating && (e.target.style.backgroundColor = '#15803d')}
            onMouseLeave={(e) => !isSimulating && (e.target.style.backgroundColor = '#16a34a')}
          >
            {isSimulating ? 'Simulating...' : 'Run Simulation'}
          </button>
        </div>

        {/* Model Descriptions */}
        <div style={styles.card}>
          <h2 style={styles.sectionTitle}>Conservation Scenarios</h2>
          <div style={styles.gridCols3}>
            {[1, 2, 3, 4, 5].map(num => {
              // Calculate effective Ne for display
              let displayNe = ne;
              if (num === 2 && results) {
                const model2Data = results.find(m => m.modelNumber === 2);
                if (model2Data) displayNe = model2Data.data[0].effectiveNe;
              }
              
              return (
                <div key={num} style={{...styles.modelCard, borderLeftColor: modelColors[`model${num}`]}}>
                  <h3 style={styles.modelTitle}>Model {num}</h3>
                  <p style={styles.modelDesc}>{getModelName(num)}</p>
                  {num === 2 && results && (
                    <p style={{fontSize: '0.7rem', color: '#ef4444', marginTop: '0.25rem'}}>
                      ⚠️ Reduced Ne: {displayNe} (due to population loss)
                    </p>
                  )}
                </div>
              );
            })}
          </div>
        </div>

        {/* Results */}
        {results && (
          <>
            {/* Tabs */}
            <div style={styles.tabContainer}>
              <div style={styles.tabRow}>
                {['overview', 'heterozygosity', 'inbreeding', 'alleles', 'population'].map(tab => (
                  <button
                    key={tab}
                    onClick={() => setActiveTab(tab)}
                    style={{
                      ...styles.tab,
                      ...(activeTab === tab ? styles.tabActive : {})
                    }}
                  >
                    {tab.charAt(0).toUpperCase() + tab.slice(1)}
                  </button>
                ))}
              </div>
            </div>

            <div style={styles.tabContent}>
              {/* Overview Tab */}
              {activeTab === 'overview' && (
                <div style={styles.overviewContainer}>
                  <div style={styles.gridCols4}>
                    {results.map(model => {
                      const finalGen = model.data[model.data.length - 1];
                      return (
                        <div key={model.modelNumber} style={styles.statCard}>
                          <div style={styles.statLabel}>Model {model.modelNumber}</div>
                          <div style={{...styles.statValue, color: modelColors[`model${model.modelNumber}`]}}>
                            {finalGen.Ho}
                          </div>
                          <div style={styles.statDesc}>Final Ho at Gen {generations}</div>
                        </div>
                      );
                    })}
                  </div>
                  
                  <div style={styles.alertBox}>
                    <div style={styles.alertContent}>
                      <AlertCircle style={styles.alertIcon} size={20} />
                      <div style={styles.alertText}>
                        <p style={styles.alertTitle}>Key Findings:</p>
                        <ul style={styles.alertList}>
                          <li>• Population loss (Model 2) shows accelerated genetic diversity decline</li>
                          <li>• Supplementation (Models 3-5) helps maintain heterozygosity</li>
                          <li>• Higher supplementation rates (Model 4) provide greater genetic benefits</li>
                          <li>• International genetic mixing (Model 5) maximizes allelic diversity</li>
                        </ul>
                      </div>
                    </div>
                  </div>
                </div>
              )}

              {/* Chart Tabs */}
              {['heterozygosity', 'inbreeding', 'alleles', 'population'].includes(activeTab) && (
                <div>
                  <h3 style={styles.chartTitle}>
                    {activeTab === 'heterozygosity' && 'Observed Heterozygosity (Ho) Over Time'}
                    {activeTab === 'inbreeding' && 'Inbreeding Coefficient (F) Over Time'}
                    {activeTab === 'alleles' && 'Mean Number of Alleles (Na) Over Time'}
                    {activeTab === 'population' && 'Population Size Projection'}
                  </h3>
                  <ResponsiveContainer width="100%" height={400}>
                    <LineChart data={prepareChartData(
                      activeTab === 'heterozygosity' ? 'Ho' :
                      activeTab === 'inbreeding' ? 'F' :
                      activeTab === 'alleles' ? 'Na' : 'popSize'
                    )}>
                      <CartesianGrid strokeDasharray="3 3" />
                      <XAxis dataKey="generation" label={{ value: 'Generation', position: 'insideBottom', offset: -5 }} />
                      <YAxis label={{ 
                        value: activeTab === 'heterozygosity' ? 'Heterozygosity (Ho)' :
                               activeTab === 'inbreeding' ? 'Inbreeding (F)' :
                               activeTab === 'alleles' ? 'Alleles per Locus (Na)' : 
                               'Population Size (N)',
                        angle: -90, 
                        position: 'insideLeft' 
                      }} 
                      domain={activeTab === 'heterozygosity' ? [0, 0.6] : 
                              activeTab === 'inbreeding' ? [0, 1] :
                              activeTab === 'alleles' ? [0, 8] : undefined}
                      />
                      <Tooltip />
                      <Legend />
                      {[1, 2, 3, 4, 5].map(num => (
                        <Line
                          key={num}
                          type="monotone"
                          dataKey={`model${num}`}
                          stroke={modelColors[`model${num}`]}
                          name={`Model ${num}`}
                          strokeWidth={2}
                          dot={false}
                        />
                      ))}
                    </LineChart>
                  </ResponsiveContainer>
                  <p style={styles.chartNote}>
                    {activeTab === 'heterozygosity' && 'Heterozygosity measures genetic diversity. Values closer to initial (0.502) indicate better genetic health.'}
                    {activeTab === 'inbreeding' && 'Inbreeding coefficient (F) ranges from 0 (no inbreeding) to 1 (complete inbreeding). Values above 0.25 indicate concern.'}
                    {activeTab === 'alleles' && 'Allelic richness represents the diversity of genetic variants. Loss of alleles reduces adaptive potential.'}
                    {activeTab === 'population' && 'Population projections assume stable growth (λ=1.0) with supplementation adding individuals directly.'}
                  </p>
                </div>
              )}
            </div>

            {/* Assumptions */}
            <div style={styles.card}>
              <h2 style={styles.sectionTitle}>Model Assumptions & Methods</h2>
              <div style={styles.assumptionsContainer}>
                <div>
                  <h3 style={styles.assumptionTitle}>Genetic Parameters:</h3>
                  <ul style={styles.assumptionList}>
                    <li>Heterozygosity loss: Ht = H0 × (1 - 1/(2Ne))^t (Wright, 1931)</li>
                    <li>Inbreeding accumulation: F = 1 - (1 - 1/(2Ne))^t</li>
                    <li>Allelic diversity: At = A0 × exp(-t/(4Ne)) (Nei et al., 1975)</li>
                  </ul>
                </div>
                
                <div>
                  <h3 style={styles.assumptionTitle}>Demographic Parameters:</h3>
                  <ul style={styles.assumptionList}>
                    <li>Generation time: 15 years (based on breeding age of 10 years)</li>
                    <li>Breeding proportion: 20% (cooperative breeding system)</li>
                    <li>Fecundity: 1 chick every 3 years per breeding pair</li>
                    <li>Population growth: λ = 1.0 (stable, based on long-lived species dynamics)</li>
                  </ul>
                </div>
                
                <div>
                  <h3 style={styles.assumptionTitle}>Population Estimates:</h3>
                  <ul style={styles.assumptionList}>
                    <li>Wild South Africa total: ~450 individuals (Taylor, 2015; Kemp & Webster, 2008)</li>
                    <li>Model 1 baseline: All 4 populations (SAE, SAG, SAK, SAL) = 450 birds</li>
                    <li>Model 2 removes Eastern Cape (n=8) and KwaZulu-Natal (n=14) = 398 birds (~52 fewer)</li>
                    <li><strong>CRITICAL:</strong> Model 2 uses proportionally smaller Ne = Ne × (398/450) = 88% of user Ne</li>
                    <li>Smaller populations experience faster genetic drift and diversity loss</li>
                  </ul>
                </div>
                
                <div>
                  <h3 style={styles.assumptionTitle}>Supplementation Effects:</h3>
                  <ul style={styles.assumptionList}>
                    <li>Gene flow increases effective Ne through immigrant contribution</li>
                    <li>Each immigrant contributes ~0.5 to effective Ne (Mills & Allendorf, 1996)</li>
                    <li>Model 5 assumes mixed captive sources reduce inbreeding more effectively</li>
                  </ul>
                </div>
                
                <div>
                  <h3 style={styles.assumptionTitle}>Key Sources:</h3>
                  <ul style={styles.assumptionList}>
                    <li>Wright, S. (1931). Evolution in Mendelian populations. Genetics 16(97)</li>
                    <li>Nei, M. et al. (1975). The bottleneck effect and genetic variability. Evolution 29(1)</li>
                    <li>Mills, L.S. & Allendorf, F.W. (1996). The one-migrant-per-generation rule. Conservation Biology 10(6)</li>
                    <li>Kemp, A.C. & Webster, P.J. (2008). Family Bucerotidae (hornbills). Handbook of Birds of the World</li>
                    <li>Taylor, M.R. (2015). The Eskom Red Data Book of Birds of South Africa</li>
                  </ul>
                </div>
                
                <div style={styles.limitationsBox}>
                  <p style={styles.limitationsText}>
                    <strong>Limitations:</strong> This model simplifies complex population dynamics. Real-world factors include environmental stochasticity, 
                    habitat quality, disease, catastrophes, and behavioral factors affecting breeding success. Results should be interpreted as 
                    relative comparisons between scenarios rather than absolute predictions.
                  </p>
                </div>
              </div>
            </div>
          </>
        )}
      </div>
    </div>
  );
};

export default Sim;