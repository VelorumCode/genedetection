from flask import Flask, render_template, request, jsonify
from main import GeneDatabase, DiseaseAnalyzer # Import from your refactored main.py

app = Flask(__name__)

# --- Initialization ---
# Create global instances to load data only once
try:
    gene_db = GeneDatabase() # Loads data from database.json on init
    analyzer = DiseaseAnalyzer(gene_db)
    if not gene_db.diseases_data:
      app.logger.warning("Disease database is empty or failed to load. Analysis may not work.")

except Exception as e:
    app.logger.error(f"Failed to initialize GeneDatabase or DiseaseAnalyzer: {e}")
    # Depending on severity, you might want to exit or prevent the app from running
    gene_db = None
    analyzer = None

# --- Routes ---
@app.route('/')
def index():
    """Serve the main HTML page"""
    # Check if initialization failed
    if not gene_db or not analyzer:
         return "Error: Application failed to initialize. Please check logs.", 500
    return render_template('index.html')

@app.route('/analyze', methods=['POST'])
def analyze_dna():
    """Handle DNA analysis requests"""
    # Check if initialization failed
    if not gene_db or not analyzer:
         return jsonify({'error': 'Application not initialized properly'}), 500

    try:
        # --- Input Validation ---
        dna_sequence = request.form.get('dna_sequence', '').strip()
        if not dna_sequence:
            return jsonify({'error': 'DNA sequence is required'}), 400

        # Basic validation (allow only ATCG, ignore case for now)
        dna_sequence = dna_sequence.upper()
        if not all(base in 'ATCG' for base in dna_sequence):
             # More robust approach might allow 'N' or other chars, but for this
             # example we require pure ATCG after converting to upper case.
             return jsonify({'error': 'Invalid characters in DNA sequence. Use only A, T, C, G'}), 400

        if len(dna_sequence) < 10: # Arbitrary minimum length
             return jsonify({'error': 'DNA sequence seems too short for meaningful analysis'}), 400

        # --- Optional Parameters ---
        age_str = request.form.get('age', '').strip()
        age: Optional[int] = None
        if age_str:
            try:
                age = int(age_str)
                if not (0 <= age <= 120): # Realistic age range
                    return jsonify({'error': 'Age must be between 0 and 120'}), 400
            except ValueError:
                return jsonify({'error': 'Age must be a valid whole number'}), 400

        gender = request.form.get('gender', '').strip().upper()
        if gender and gender not in ['M', 'F']:
            return jsonify({'error': 'Gender must be M or F (or leave blank)'}), 400
        if not gender:
            gender = None # Ensure gender is None if blank

        # --- Perform Analysis ---
        app.logger.info(f"Analyzing sequence (length {len(dna_sequence)}) with age={age}, gender={gender}")
        matches = gene_db.search_markers(dna_sequence)

        if not matches:
            app.logger.info("No matching disease markers found.")
            # Return success, but indicate no findings
            return jsonify({
                'message': 'No known disease markers found in the provided sequence.',
                'results': {}
                }), 200

        app.logger.info(f"Found potential markers for diseases: {list(matches.keys())}")
        analysis_results = analyzer.calculate_risk_score(matches, patient_age=age, gender=gender)

        # --- Return Results ---
        return jsonify({
             'message': 'Analysis complete. See results below.',
             'results': analysis_results
             })

    except Exception as e:
        # Log the full error for debugging
        app.logger.error(f"An unexpected error occurred during analysis: {str(e)}", exc_info=True)
        # Return a generic error to the user
        return jsonify({'error': 'An internal error occurred during analysis. Please try again later.'}), 500

# --- Run Application ---
if __name__ == '__main__':
    # Set debug=False for production
    app.run(debug=True)