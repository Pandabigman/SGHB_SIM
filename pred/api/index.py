from flask import Flask, render_template, request, jsonify
import os
import sys
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Add parent directory to path for package imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from genetics import LOCI, load_and_process_genetic_data
from models import run_model

# Get base directory for path resolution
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

app = Flask(
    __name__,
    template_folder=os.path.join(BASE_DIR, "templates"),
    static_folder=os.path.join(BASE_DIR, "static"),
)


# LOAD GENETIC DATA FROM CSVs


# Paths to CSV files
WILD_CSV = os.path.join(BASE_DIR, "data", "blue.csv")
CAPTIVE_CSV = os.path.join(BASE_DIR, "data", "red.csv")

# Global variable to store processed genetic data
GENETIC_DATA = None


def initialize_genetic_data():
    """Load and process CSV data on startup"""
    global GENETIC_DATA

    if os.path.exists(WILD_CSV) and os.path.exists(CAPTIVE_CSV):
        logger.info("Loading genetic data from CSVs...")
        GENETIC_DATA = load_and_process_genetic_data(WILD_CSV, CAPTIVE_CSV)
        logger.info(f"Loaded {len(GENETIC_DATA['summaries'])} populations")
        logger.info(
            f"Wild population Ho: {GENETIC_DATA['summaries']['wild_all']['Ho']:.4f}"
        )
    else:
        logger.warning("CSV files not found. Using default values.")
        GENETIC_DATA = None


# API ROUTES


@app.route("/")
def index():
    """Serve the main page"""
    return render_template("index.html")


@app.route("/api/simulate", methods=["POST"])
def simulate():
    """Run simulation for specified model"""
    try:
        data = request.get_json()

        Ne = int(data.get("Ne", 500))
        generations = int(data.get("generations", 50))
        lambda_val = float(data.get("lambda", 1.0))
        model = int(data.get("model", 1))
        stochastic = bool(data.get("stochastic", False))

        # Validate inputs
        if not 375 <= Ne <= 625:
            return jsonify({"error": "Ne must be between 375 and 625"}), 400
        if not 10 <= generations <= 100:
            return jsonify({"error": "Generations must be between 10 and 100"}), 400
        if not 0.5 <= lambda_val <= 1.5:
            return jsonify({"error": "Lambda must be between 0.5 and 1.5"}), 400
        if not 1 <= model <= 6:
            return jsonify({"error": "Model must be between 1 and 6"}), 400

        # Run the model using the models package
        results = run_model(
            model_num=model,
            Ne=Ne,
            generations=generations,
            lambda_val=lambda_val,
            stochastic=stochastic,
            genetic_data=GENETIC_DATA
        )

        return jsonify(results)

    except Exception as e:
        import traceback
        traceback.print_exc()
        logger.error(f"Simulation error: {str(e)}")
        return jsonify({"error": str(e)}), 500


@app.route("/api/data/info", methods=["GET"])
def get_data_info():
    """Get information about genetic data"""
    if GENETIC_DATA:
        return jsonify(
            {
                "data_source": "CSV",
                "populations": {
                    name: {
                        "sample_size": summary["sample_size"],
                        "Ho": summary["Ho"],
                        "He": summary["He"],
                        "Na": summary["Na"],
                        "FIS": summary["FIS"],
                    }
                    for name, summary in GENETIC_DATA["summaries"].items()
                },
                "lost_alleles_count": sum(
                    len(alleles)
                    for alleles in GENETIC_DATA["lost_alleles_ec_kzn"].values()
                ),
                "novel_alleles": {
                    "paaza": sum(
                        len(a) for a in GENETIC_DATA["novel_alleles"]["paaza"].values()
                    ),
                    "aza": sum(
                        len(a) for a in GENETIC_DATA["novel_alleles"]["aza"].values()
                    ),
                    "eaza": sum(
                        len(a) for a in GENETIC_DATA["novel_alleles"]["eaza"].values()
                    ),
                },
            }
        )
    else:
        return jsonify(
            {
                "data_source": "default",
                "message": "Using default genetic parameters. Upload CSVs for real data.",
            }
        )


# Initialize data on import
initialize_genetic_data()

# For Vercel serverless functions
app = app

if __name__ == "__main__":
    app.run(debug=True, port=5001)
