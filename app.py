from flask import Flask, render_template, request, jsonify
import random
import math
from main import FusionTree, DiseaseAnalyzer, initialize_disease_database

app = Flask(__name__)

# Create global instances to avoid reinitializing for each request
fusion_tree = FusionTree()
analyzer = DiseaseAnalyzer(fusion_tree)

# Initialize the disease database during startup
diseases = initialize_disease_database()
for disease, info in diseases.items():
    fusion_tree.add_gene(
        disease=disease,
        gene_sequences=info["sequences"],
        prevalence=info["prevalence"],
        age_risk=info["age_risk"],
        gender_risk=info["gender_risk"]
    )

@app.route('/')
def index():
    """Serve the main HTML page"""
    return render_template('index.html')

@app.route('/analyze', methods=['POST'])
def analyze_dna():
    """Handle DNA analysis requests"""
    try:
        # Get and validate input data
        dna_sequence = request.form.get('dna_sequence', '').strip().upper()
        if not dna_sequence:
            return jsonify({'error': 'DNA sequence is required'}), 400

        # Validate DNA sequence (should only contain A, T, C, G)
        if not all(base in 'ATCG' for base in dna_sequence):
            return jsonify({'error': 'Invalid DNA sequence. Use only A, T, C, G'}), 400

        # Get optional parameters
        try:
            age = int(request.form.get('age', '')) if request.form.get('age') else None
        except ValueError:
            return jsonify({'error': 'Age must be a valid number'}), 400

        gender = request.form.get('gender', '').upper()
        if gender and gender not in ['M', 'F']:
            return jsonify({'error': 'Gender must be M or F'}), 400

        # Perform analysis
        matches = fusion_tree.search_gene(dna_sequence)
        if not matches:
            return jsonify({'message': 'No significant genetic matches found'}), 200

        analysis_results = analyzer.calculate_probability(matches, patient_age=age, gender=gender)
        
        return jsonify(analysis_results)

    except Exception as e:
        app.logger.error(f"Analysis error: {str(e)}")
        return jsonify({'error': 'An error occurred during analysis'}), 500

if __name__ == '__main__':
    app.run(debug=True)