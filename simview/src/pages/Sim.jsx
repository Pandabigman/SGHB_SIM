import React, { useState, useEffect } from "react";
import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  Legend,
  ResponsiveContainer,
} from "recharts";
import { Info, AlertCircle } from "lucide-react";
import { runSimulationModels, getModelName } from "../js/simulation";

const Sim = () => {
  // State for user inputs
  const [ne, setNe] = useState(100);
  const [generations, setGenerations] = useState(50);

  // State for simulation results and UI
  const [isSimulating, setIsSimulating] = useState(false);
  const [results, setResults] = useState(null);
  const [activeTab, setActiveTab] = useState("overview");

  const runSimulation = () => {
    setIsSimulating(true);
    // The complex calculations are now in an imported function
    const models = runSimulationModels({ ne, generations });
    setResults(models);
    setIsSimulating(false);
  };

  // Prepare data for combined charts
  const prepareChartData = (metric) => {
    if (!results) return [];

    const combinedData = [];
    const maxGen = generations;

    for (let gen = 0; gen <= maxGen; gen++) {
      const dataPoint = { generation: gen };

      results.forEach((model) => {
        const genData = model.data.find((d) => d.generation === gen);
        if (genData) {
          dataPoint[`model${model.modelNumber}`] = parseFloat(genData[metric]);
        }
      });

      combinedData.push(dataPoint);
    }

    return combinedData;
  };

  const modelColors = {
    model1: "#2563eb", // Blue
    model2: "#dc2626", // Red
    model3: "#16a34a", // Green
    model4: "#9333ea", // Purple
    model5: "#ea580c", // Orange
  };

  useEffect(() => {
    runSimulation();
  }, []);

  return (
    <div className="min-h-screen bg-gradient-to-br from-green-50 to-blue-50 p-6">
      <div className="max-w-7xl mx-auto">
        {/* Header */}
        <div className="bg-white rounded-lg shadow-lg p-6 mb-6">
          <div className="flex items-start gap-4">
            <div className="flex-1">
              <h1 className="text-3xl font-bold text-gray-800 mb-2">
                Southern Ground Hornbill Population Genetics Simulator
              </h1>
              <p className="text-gray-600 mb-4">
                Predictive modeling of genetic diversity and population
                viability over 50 generations (750 years)
              </p>
              <div className="flex items-center gap-2 text-sm text-gray-500">
                <Info size={16} />
                <span>
                  Based on microsatellite data from 199 wild and 145 captive
                  individuals
                </span>
              </div>
            </div>
          </div>
        </div>

        {/* Controls */}
        <div className="bg-white rounded-lg shadow-lg p-6 mb-6">
          <h2 className="text-xl font-semibold mb-4">Simulation Parameters</h2>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            <div>
              <label className="block text-sm font-medium text-gray-700 mb-2">
                Effective Population Size (Ne): {ne}
              </label>
              <input
                type="range"
                min="20"
                max="200"
                value={ne}
                onChange={(e) => setNe(parseInt(e.target.value))}
                className="w-full h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer"
              />
              <div className="flex justify-between text-xs text-gray-500 mt-1">
                <span>20</span>
                <span>100</span>
                <span>200</span>
              </div>
            </div>

            <div>
              <label className="block text-sm font-medium text-gray-700 mb-2">
                Generations to Simulate: {generations}
              </label>
              <input
                type="range"
                min="10"
                max="100"
                step="10"
                value={generations}
                onChange={(e) => setGenerations(parseInt(e.target.value))}
                className="w-full h-2 bg-gray-200 rounded-lg appearance-none cursor-pointer"
              />
              <div className="flex justify-between text-xs text-gray-500 mt-1">
                <span>10 (150 yrs)</span>
                <span>50 (750 yrs)</span>
                <span>100 (1500 yrs)</span>
              </div>
            </div>
          </div>

          <button
            onClick={runSimulation}
            disabled={isSimulating}
            className="mt-6 w-full md:w-auto px-6 py-3 bg-green-600 hover:bg-green-700 text-white font-medium rounded-lg transition-colors disabled:bg-gray-400"
          >
            {isSimulating ? "Simulating..." : "Run Simulation"}
          </button>
        </div>

        {/* Model Descriptions */}
        <div className="bg-white rounded-lg shadow-lg p-6 mb-6">
          <h2 className="text-xl font-semibold mb-4">Conservation Scenarios</h2>
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
            {[1, 2, 3, 4, 5].map((num) => (
              <div
                key={num}
                className="border-l-4 pl-4"
                style={{ borderColor: modelColors[`model${num}`] }}
              >
                <h3 className="font-semibold text-sm mb-1">Model {num}</h3>
                <p className="text-xs text-gray-600">{getModelName(num)}</p>
              </div>
            ))}
          </div>
        </div>

        {/* Results */}
        {results && (
          <>
            {/* Tabs */}
            <div className="bg-white rounded-t-lg shadow-lg">
              <div className="flex border-b overflow-x-auto">
                {[
                  "overview",
                  "heterozygosity",
                  "inbreeding",
                  "alleles",
                  "population",
                ].map((tab) => (
                  <button
                    key={tab}
                    onClick={() => setActiveTab(tab)}
                    className={`px-6 py-3 font-medium transition-colors ${
                      activeTab === tab
                        ? "text-green-600 border-b-2 border-green-600"
                        : "text-gray-600 hover:text-gray-800"
                    }`}
                  >
                    {tab.charAt(0).toUpperCase() + tab.slice(1)}
                  </button>
                ))}
              </div>
            </div>

            <div className="bg-white rounded-b-lg shadow-lg p-6">
              {/* Overview Tab */}
              {activeTab === "overview" && (
                <div className="space-y-6">
                  <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                    {results.map((model) => {
                      const finalGen = model.data[model.data.length - 1];
                      return (
                        <div
                          key={model.modelNumber}
                          className="border rounded-lg p-4"
                        >
                          <div className="text-xs font-medium text-gray-500 mb-1">
                            Model {model.modelNumber}
                          </div>
                          <div
                            className="text-2xl font-bold mb-2"
                            style={{
                              color: modelColors[`model${model.modelNumber}`],
                            }}
                          >
                            {finalGen.Ho}
                          </div>
                          <div className="text-xs text-gray-600">
                            Final Ho at Gen {generations}
                          </div>
                        </div>
                      );
                    })}
                  </div>

                  <div className="bg-blue-50 border border-blue-200 rounded-lg p-4">
                    <div className="flex items-start gap-2">
                      <AlertCircle
                        size={20}
                        className="text-blue-600 mt-0.5 flex-shrink-0"
                      />
                      <div className="text-sm text-blue-900">
                        <p className="font-semibold mb-1">Key Findings:</p>
                        <ul className="space-y-1 text-xs">
                          <li>
                            • Population loss (Model 2) shows accelerated
                            genetic diversity decline
                          </li>
                          <li>
                            • Supplementation (Models 3-5) helps maintain
                            heterozygosity
                          </li>
                          <li>
                            • Higher supplementation rates (Model 4) provide
                            greater genetic benefits
                          </li>
                          <li>
                            • International genetic mixing (Model 5) maximizes
                            allelic diversity
                          </li>
                        </ul>
                      </div>
                    </div>
                  </div>
                </div>
              )}

              {/* Heterozygosity Tab */}
              {activeTab === "heterozygosity" && (
                <div>
                  <h3 className="text-lg font-semibold mb-4">
                    Observed Heterozygosity (Ho) Over Time
                  </h3>
                  <ResponsiveContainer width="100%" height={400}>
                    <LineChart data={prepareChartData("Ho")}>
                      <CartesianGrid strokeDasharray="3 3" />
                      <XAxis
                        dataKey="generation"
                        label={{
                          value: "Generation",
                          position: "insideBottom",
                          offset: -5,
                        }}
                      />
                      <YAxis
                        label={{
                          value: "Heterozygosity (Ho)",
                          angle: -90,
                          position: "insideLeft",
                        }}
                        domain={[0, 0.6]}
                      />
                      <Tooltip />
                      <Legend />
                      {[1, 2, 3, 4, 5].map((num) => (
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
                  <p className="text-sm text-gray-600 mt-4">
                    Heterozygosity measures genetic diversity. Values closer to
                    initial (0.502) indicate better genetic health.
                  </p>
                </div>
              )}

              {/* Inbreeding Tab */}
              {activeTab === "inbreeding" && (
                <div>
                  <h3 className="text-lg font-semibold mb-4">
                    Inbreeding Coefficient (F) Over Time
                  </h3>
                  <ResponsiveContainer width="100%" height={400}>
                    <LineChart data={prepareChartData("F")}>
                      <CartesianGrid strokeDasharray="3 3" />
                      <XAxis
                        dataKey="generation"
                        label={{
                          value: "Generation",
                          position: "insideBottom",
                          offset: -5,
                        }}
                      />
                      <YAxis
                        label={{
                          value: "Inbreeding (F)",
                          angle: -90,
                          position: "insideLeft",
                        }}
                        domain={[0, 1]}
                      />
                      <Tooltip />
                      <Legend />
                      {[1, 2, 3, 4, 5].map((num) => (
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
                  <p className="text-sm text-gray-600 mt-4">
                    Inbreeding coefficient (F) ranges from 0 (no inbreeding) to
                    1 (complete inbreeding). Values above 0.25 indicate concern.
                  </p>
                </div>
              )}

              {/* Alleles Tab */}
              {activeTab === "alleles" && (
                <div>
                  <h3 className="text-lg font-semibold mb-4">
                    Mean Number of Alleles (Na) Over Time
                  </h3>
                  <ResponsiveContainer width="100%" height={400}>
                    <LineChart data={prepareChartData("Na")}>
                      <CartesianGrid strokeDasharray="3 3" />
                      <XAxis
                        dataKey="generation"
                        label={{
                          value: "Generation",
                          position: "insideBottom",
                          offset: -5,
                        }}
                      />
                      <YAxis
                        label={{
                          value: "Alleles per Locus (Na)",
                          angle: -90,
                          position: "insideLeft",
                        }}
                        domain={[0, 8]}
                      />
                      <Tooltip />
                      <Legend />
                      {[1, 2, 3, 4, 5].map((num) => (
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
                  <p className="text-sm text-gray-600 mt-4">
                    Allelic richness represents the diversity of genetic
                    variants. Loss of alleles reduces adaptive potential.
                  </p>
                </div>
              )}

              {/* Population Tab */}
              {activeTab === "population" && (
                <div>
                  <h3 className="text-lg font-semibold mb-4">
                    Population Size Projection
                  </h3>
                  <ResponsiveContainer width="100%" height={400}>
                    <LineChart data={prepareChartData("popSize")}>
                      <CartesianGrid strokeDasharray="3 3" />
                      <XAxis
                        dataKey="generation"
                        label={{
                          value: "Generation",
                          position: "insideBottom",
                          offset: -5,
                        }}
                      />
                      <YAxis
                        label={{
                          value: "Population Size (N)",
                          angle: -90,
                          position: "insideLeft",
                        }}
                      />
                      <Tooltip />
                      <Legend />
                      {[1, 2, 3, 4, 5].map((num) => (
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
                  <p className="text-sm text-gray-600 mt-4">
                    Population projections assume stable growth (λ=1.0) with
                    supplementation adding individuals directly.
                  </p>
                </div>
              )}
            </div>

            {/* Assumptions */}
            <div className="bg-white rounded-lg shadow-lg p-6 mt-6">
              <h2 className="text-xl font-semibold mb-4">
                Model Assumptions & Methods
              </h2>
              <div className="space-y-4 text-sm text-gray-700">
                <div>
                  <h3 className="font-semibold mb-2">Genetic Parameters:</h3>
                  <ul className="space-y-1 ml-4 list-disc">
                    <li>
                      Heterozygosity loss: Ht = H0 × (1 - 1/(2Ne))^t (Wright,
                      1931)
                    </li>
                    <li>Inbreeding accumulation: F = 1 - (1 - 1/(2Ne))^t</li>
                    <li>
                      Allelic diversity: At = A0 × exp(-t/(4Ne)) (Nei et al.,
                      1975)
                    </li>
                  </ul>
                </div>

                <div>
                  <h3 className="font-semibold mb-2">
                    Demographic Parameters:
                  </h3>
                  <ul className="space-y-1 ml-4 list-disc">
                    <li>
                      Generation time: 15 years (based on breeding age of 10
                      years)
                    </li>
                    <li>
                      Breeding proportion: 20% (cooperative breeding system)
                    </li>
                    <li>Fecundity: 1 chick every 3 years per breeding pair</li>
                    <li>
                      Population growth: λ = 1.0 (stable, based on long-lived
                      species dynamics)
                    </li>
                  </ul>
                </div>

                <div>
                  <h3 className="font-semibold mb-2">Population Estimates:</h3>
                  <ul className="space-y-1 ml-4 list-disc">
                    <li>
                      Wild South Africa total: ~450 individuals (Taylor, 2015;
                      Kemp & Webster, 2008)
                    </li>
                    <li>
                      Proportional allocation based on sample sizes from study
                      regions
                    </li>
                    <li>
                      Model 2 removes Eastern Cape (n=8) and KwaZulu-Natal
                      (n=14) contributions
                    </li>
                  </ul>
                </div>

                <div>
                  <h3 className="font-semibold mb-2">
                    Supplementation Effects:
                  </h3>
                  <ul className="space-y-1 ml-4 list-disc">
                    <li>
                      Gene flow increases effective Ne through immigrant
                      contribution
                    </li>
                    <li>
                      Each immigrant contributes ~0.5 to effective Ne (Mills &
                      Allendorf, 1996)
                    </li>
                    <li>
                      Model 5 assumes mixed captive sources reduce inbreeding
                      more effectively
                    </li>
                  </ul>
                </div>

                <div>
                  <h3 className="font-semibold mb-2">Key Sources:</h3>
                  <ul className="space-y-1 ml-4 list-disc">
                    <li>
                      Wright, S. (1931). Evolution in Mendelian populations.
                      Genetics 16(97)
                    </li>
                    <li>
                      Nei, M. et al. (1975). The bottleneck effect and genetic
                      variability. Evolution 29(1)
                    </li>
                    <li>
                      Mills, L.S. & Allendorf, F.W. (1996). The
                      one-migrant-per-generation rule. Conservation Biology
                      10(6)
                    </li>
                    <li>
                      Kemp, A.C. & Webster, P.J. (2008). Family Bucerotidae
                      (hornbills). Handbook of Birds of the World
                    </li>
                    <li>
                      Taylor, M.R. (2015). The Eskom Red Data Book of Birds of
                      South Africa
                    </li>
                  </ul>
                </div>

                <div className="bg-yellow-50 border border-yellow-200 rounded-lg p-3">
                  <p className="text-xs">
                    <strong>Limitations:</strong> This model simplifies complex
                    population dynamics. Real-world factors include
                    environmental stochasticity, habitat quality, disease,
                    catastrophes, and behavioral factors affecting breeding
                    success. Results should be interpreted as relative
                    comparisons between scenarios rather than absolute
                    predictions.
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
