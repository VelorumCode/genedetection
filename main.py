import json
from typing import Dict, List, Tuple, Optional, Any

# --- Constants ---
DATABASE_FILE = 'database.json'
MARKER_MATCH_THRESHOLD = 0.8 # Example: Require 80% match if using alignment, or exact match if using 'in'

# --- Helper Functions ---

def load_disease_database(filepath: str = DATABASE_FILE) -> Dict[str, Dict[str, Any]]:
    """Loads the disease database from a JSON file."""
    try:
        with open(filepath, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        print(f"Error: Database file '{filepath}' not found.")
        return {}
    except json.JSONDecodeError:
        print(f"Error: Could not decode JSON from '{filepath}'.")
        return {}

# --- Core Classes ---

class GeneDatabase:
    """
    Stores disease information, associated genetic markers, and metadata.
    Provides methods to search for markers within a given DNA sequence.
    """
    def __init__(self, database_path: str = DATABASE_FILE):
        self.diseases_data = load_disease_database(database_path)
        if not self.diseases_data:
            print("Warning: Gene database is empty or failed to load.")

    def get_disease_data(self, disease: str) -> Optional[Dict[str, Any]]:
        """Retrieve metadata for a specific disease."""
        return self.diseases_data.get(disease)

    def search_markers(self, input_sequence: str) -> Dict[str, List[str]]:
        """
        Searches for known disease markers within the input DNA sequence.

        Args:
            input_sequence: The DNA sequence provided by the user.

        Returns:
            A dictionary where keys are disease names and values are lists
            of marker sequences found in the input sequence for that disease.
            Returns an empty dictionary if no matches are found.
        """
        matches: Dict[str, List[str]] = {}
        if not self.diseases_data or not input_sequence:
            return matches # Return empty if no data or no sequence

        # Ensure input is valid DNA
        input_sequence = input_sequence.upper()
        if not all(base in 'ATCG' for base in input_sequence):
             # In a real scenario, might handle 'N' or other IUPAC codes
             print("Warning: Input sequence contains non-ATCG characters.")
             # Decide whether to filter them out or reject the sequence
             input_sequence = ''.join(filter(lambda base: base in 'ATCG', input_sequence))
             if not input_sequence: return matches # Return empty if only invalid chars

        for disease, info in self.diseases_data.items():
            found_markers_for_disease = []
            markers = info.get("markers", [])
            for marker in markers:
                marker = marker.upper() # Ensure marker is uppercase
                # --- Simple Substring Matching ---
                # Check if the marker sequence exists within the input sequence
                if marker in input_sequence:
                    found_markers_for_disease.append(marker)

                # --- Alternative: Basic Similarity (Less Recommended for Variants) ---
                # If you wanted to keep a similarity score, you'd need a function:
                # similarity = calculate_similarity(input_sequence, marker) # Needs a robust function
                # if similarity >= MARKER_MATCH_THRESHOLD:
                #    found_markers_for_disease.append((marker, similarity)) # Store marker and score

            if found_markers_for_disease:
                matches[disease] = found_markers_for_disease

        return matches

class DiseaseAnalyzer:
    """
    Analyzes the results from GeneDatabase to calculate risk scores
    based on found markers and patient demographic data.
    """
    def __init__(self, gene_database: GeneDatabase):
        self.db = gene_database

    def calculate_risk_score(self,
                             matches: Dict[str, List[str]],
                             patient_age: Optional[int] = None,
                             gender: Optional[str] = None) -> Dict[str, Dict[str, Any]]:
        """
        Calculates a risk score for each matched disease based on genetic
        matches, disease prevalence, and demographic factors.

        Args:
            matches: Dictionary of diseases and the markers found for them.
            patient_age: Optional age of the patient.
            gender: Optional gender of the patient ('M' or 'F').

        Returns:
            A dictionary where keys are disease names and values contain
            details about the risk calculation (markers found, risk factors,
            and a final 'risk_score').
        """
        if not matches:
            return {}

        analysis_results: Dict[str, Dict[str, Any]] = {}
        gender = gender.upper() if gender else None

        for disease, found_markers in matches.items():
            metadata = self.db.get_disease_data(disease)
            if not metadata:
                continue # Skip if no metadata found for this disease

            # --- Base Risk Factors ---
            # 1. Genetic Component: Base risk increases if *any* marker is found.
            #    (Could be more complex, e.g., weighting different markers)
            genetic_factor = 1.0 # Base factor if a marker is found

            # 2. Prevalence Factor: Use disease prevalence. Lower prevalence might mean
            #    a found marker is more significant, or could just scale the base risk.
            #    Let's use it as a scaler. Handle potential missing data.
            prevalence_factor = metadata.get('prevalence', 0.0001) # Default low prevalence

            # 3. Age Adjustment
            age_multiplier = 1.0
            if patient_age is not None:
                age_risks = metadata.get('age_risk', {})
                for age_range, multiplier in age_risks.items():
                    try:
                        start_str, end_str = age_range.split('-')
                        start = int(start_str)
                        end = int(end_str)
                        if start <= patient_age <= end:
                            age_multiplier = float(multiplier)
                            break
                    except ValueError:
                        print(f"Warning: Invalid age range format '{age_range}' for disease '{disease}'")
                        continue # Skip this invalid range

            # 4. Gender Adjustment
            gender_multiplier = 1.0
            if gender and gender in ['M', 'F']:
                gender_risks = metadata.get('gender_risk', {})
                gender_multiplier = float(gender_risks.get(gender, 1.0))

            # --- Calculate Final Risk Score ---
            # Combine factors multiplicatively (this is a simplification)
            # A higher score indicates higher relative risk based on the model.
            # NOTE: This score is NOT a probability.
            risk_score = (genetic_factor *
                          prevalence_factor *
                          age_multiplier *
                          gender_multiplier)

            # Scale score for better readability (optional)
            # Multiplying by a large number avoids tiny decimals
            display_score = risk_score * 100000

            analysis_results[disease] = {
                'markers_found': found_markers,
                'description': metadata.get('description', 'No description available.'),
                'risk_factors': {
                   'genetic_component': genetic_factor, # Simplified: 1 if marker found
                   'prevalence': prevalence_factor,
                   'age_multiplier': age_multiplier,
                   'gender_multiplier': gender_multiplier
                },
                'calculated_risk_score': display_score # This is a relative score, not probability
            }

        # Optional: Sort results by risk score (descending)
        sorted_results = dict(sorted(analysis_results.items(),
                                     key=lambda item: item[1]['calculated_risk_score'],
                                     reverse=True))

        return sorted_results