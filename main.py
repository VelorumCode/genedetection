import random
import math
from typing import Dict, List, Tuple, Optional
import numpy as np
from datetime import datetime

class FusionTree:
    def __init__(self):
        self.genes = {}
        self.disease_metadata = {}
        
    def add_gene(self, disease: str, gene_sequences: List[str], prevalence: float = None, 
                 age_risk: Dict = None, gender_risk: Dict = None):
        """
        Add gene sequences to the Fusion Tree for a specific disease with metadata.
        
        Args:
            disease: Name of the disease
            gene_sequences: List of DNA sequences associated with the disease
            prevalence: Global prevalence of the disease (percentage)
            age_risk: Dictionary of age ranges and their risk multipliers
            gender_risk: Dictionary of gender-specific risk multipliers
        """
        self.genes[disease] = gene_sequences
        
        self.disease_metadata[disease] = {
            'prevalence': prevalence or 0.01,
            'age_risk': age_risk or {'0-100': 1.0},
            'gender_risk': gender_risk or {'M': 1.0, 'F': 1.0},
            'mutation_rate': random.uniform(0.001, 0.01)
        }
    
    def search_gene(self, sequence: str) -> Dict[str, List[Tuple[str, float]]]:
        """
        Search for matching gene sequences and return matches with similarity scores.
        """
        matches = {}
        
        for disease, gene_sequences in self.genes.items():
            disease_matches = []
            
            for gene in gene_sequences:
                similarity = calculate_sequence_similarity(sequence, gene)
                if similarity > 0.6:  # Threshold for considering a match
                    disease_matches.append((gene, similarity))
            
            if disease_matches:
                matches[disease] = disease_matches
        
        return matches

# Rest of the code remains the same, including DiseaseAnalyzer class, initialize_disease_database,
# calculate_sequence_similarity, and main function...

def calculate_sequence_similarity(seq1: str, seq2: str) -> float:
    """
    Calculate the similarity between two DNA sequences using a more sophisticated algorithm.
    Returns a similarity score between 0 and 1.
    """
    if not seq1 or not seq2:
        return 0.0
    
    # Ensure sequences are of equal length for comparison
    min_length = min(len(seq1), len(seq2))
    seq1 = seq1[:min_length]
    seq2 = seq2[:min_length]
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    similarity = matches / min_length
    
    # Apply weighted scoring for consecutive matches
    consecutive_bonus = 0
    current_streak = 0
    
    for a, b in zip(seq1, seq2):
        if a == b:
            current_streak += 1
            consecutive_bonus += (current_streak * 0.1)  # Bonus increases with streak length
        else:
            current_streak = 0
    
    # Normalize bonus and combine with base similarity
    final_similarity = min(1.0, similarity + (consecutive_bonus / len(seq1)))
    return final_similarity

# Rest of the code remains exactly the same...
class DiseaseAnalyzer:
    def __init__(self, fusion_tree: FusionTree):
        self.fusion_tree = fusion_tree
        
    def calculate_probability(self, matches: Dict, patient_age: int = None, 
                            gender: str = None) -> Dict[str, Dict]:
        """
        Calculate comprehensive disease probabilities with enhanced accuracy.
        """
        if not matches:
            return {}
        
        analysis_results = {}
        total_weighted_matches = 0
        
        for disease, genes in matches.items():
            metadata = self.fusion_tree.disease_metadata[disease]
            
            # Base genetic match probability
            base_probability = sum(similarity for _, similarity in genes) / len(genes)
            
            # Prevalence adjustment
            prevalence_factor = metadata['prevalence']
            
            # Age adjustment
            age_multiplier = 1.0
            if patient_age is not None:
                age_risks = metadata['age_risk']
                for age_range, multiplier in age_risks.items():
                    start, end = map(int, age_range.split('-'))
                    if start <= patient_age <= end:
                        age_multiplier = multiplier
                        break
            
            # Gender adjustment
            gender_multiplier = 1.0
            if gender and gender.upper() in metadata['gender_risk']:
                gender_multiplier = metadata['gender_risk'][gender.upper()]
            
            # Calculate weighted probability
            weighted_probability = (base_probability * 
                                 prevalence_factor * 
                                 age_multiplier * 
                                 gender_multiplier)
            
            total_weighted_matches += weighted_probability
            
            analysis_results[disease] = {
                'genetic_match': base_probability,
                'prevalence_factor': prevalence_factor,
                'age_multiplier': age_multiplier,
                'gender_multiplier': gender_multiplier,
                'weighted_probability': weighted_probability
            }
        
        # Normalize probabilities
        if total_weighted_matches > 0:
            for disease in analysis_results:
                analysis_results[disease]['final_probability'] = (
                    analysis_results[disease]['weighted_probability'] / total_weighted_matches * 100
                )
        
        return analysis_results

def initialize_disease_database() -> Dict:
    """
    Initialize database of diseases with genetic markers and risk factors.
    """
    return {
        "Cystic Fibrosis": {
            "sequences": ["ATCGTACGATC", "GCTAGCTAGCT", "CGTATCGATCG"],
            "prevalence": 0.04,
            "age_risk": {"0-20": 1.2, "21-40": 1.0, "41-100": 0.8},
            "gender_risk": {"M": 1.0, "F": 1.0}
        },
        "Sickle Cell Anemia": {
            "sequences": ["GTACGGTACGGT", "TACGGTACGGTA", "ACGGTACGGTAT"],
            "prevalence": 0.03,
            "age_risk": {"0-30": 1.3, "31-100": 1.0},
            "gender_risk": {"M": 1.0, "F": 1.0}
        },
        "Huntington's Disease": {
            "sequences": ["TACGGTACAGTC", "ACGGTACAGTCT", "CGGTACAGTCTA"],
            "prevalence": 0.01,
            "age_risk": {"0-30": 0.8, "31-50": 1.2, "51-100": 1.5},
            "gender_risk": {"M": 1.0, "F": 1.0}
        },
        "BRCA1 Breast Cancer": {
            "sequences": ["CGTAGCTAGTAC", "GTAGCTAGTACG", "TAGCTAGTACGT"],
            "prevalence": 0.12,
            "age_risk": {"0-30": 0.5, "31-50": 1.3, "51-100": 1.8},
            "gender_risk": {"M": 0.1, "F": 1.8}
        },
        "Early-Onset Alzheimer's": {
            "sequences": ["TAGCTAGTCCGA", "AGCTAGTCCGAT", "GCTAGTCCGATA"],
            "prevalence": 0.08,
            "age_risk": {"0-40": 0.3, "41-60": 1.2, "61-100": 2.0},
            "gender_risk": {"M": 0.8, "F": 1.2}
        }
    }

def main():
    print("Advanced Genetic Disease Detection System v2.0")
    print("=" * 60)
    
    # Initialize the system
    fusion_tree = FusionTree()
    analyzer = DiseaseAnalyzer(fusion_tree)
    
    # Load disease database
    diseases = initialize_disease_database()
    
    # Add diseases to the fusion tree
    for disease, info in diseases.items():
        fusion_tree.add_gene(
            disease=disease,
            gene_sequences=info["sequences"],
            prevalence=info["prevalence"],
            age_risk=info["age_risk"],
            gender_risk=info["gender_risk"]
        )
    
    print(f"\nSystem initialized with {len(diseases)} diseases and their genetic markers.")
    
    while True:
        print("\nOptions:")
        print("1. Analyze DNA sequence")
        print("2. List available diseases")
        print("3. Exit")
        
        choice = input("\nEnter your choice (1-3): ")
        
        if choice == "3":
            break
        elif choice == "2":
            print("\nAvailable Diseases for Analysis:")
            print("-" * 40)
            for i, disease in enumerate(diseases.keys(), 1):
                print(f"{i}. {disease}")
        elif choice == "1":
            sequence = input("\nEnter the DNA sequence to analyze: ").strip().upper()
            age = input("Enter patient age (or press Enter to skip): ")
            gender = input("Enter patient gender (M/F, or press Enter to skip): ").upper()
            
            try:
                age = int(age) if age else None
            except ValueError:
                print("Invalid age entered, continuing without age adjustment...")
                age = None
            
            if gender and gender not in ['M', 'F']:
                print("Invalid gender entered, continuing without gender adjustment...")
                gender = None
            
            # Search for matches
            matches = fusion_tree.search_gene(sequence)
            
            if matches:
                print("\nGenetic Analysis Results:")
                print("=" * 60)
                
                # Calculate comprehensive probabilities
                analysis = analyzer.calculate_probability(matches, age, gender)
                
                # Sort diseases by probability
                sorted_diseases = sorted(
                    analysis.items(),
                    key=lambda x: x[1]['final_probability'],
                    reverse=True
                )
                
                # Display results
                for disease, details in sorted_diseases:
                    print(f"\n{disease}:")
                    print(f"  Final Probability: {details['final_probability']:.2f}%")
                    print("  Risk Factors:")
                    print(f"    - Genetic Match: {details['genetic_match']:.2f}")
                    print(f"    - Prevalence Factor: {details['prevalence_factor']:.3f}")
                    print(f"    - Age Risk Multiplier: {details['age_multiplier']:.2f}")
                    print(f"    - Gender Risk Multiplier: {details['gender_multiplier']:.2f}")
            else:
                print("\nNo significant genetic matches found in the database.")
        else:
            print("\nInvalid choice. Please try again.")

    print("\nThank you for using the Advanced Genetic Disease Detection System!")

def calculate_sequence_similarity(seq1: str, seq2: str) -> float:
    """
    Calculate the similarity between two DNA sequences using a more sophisticated algorithm.
    Returns a similarity score between 0 and 1.
    """
    if not seq1 or not seq2:
        return 0.0
    
    # Ensure sequences are of equal length for comparison
    min_length = min(len(seq1), len(seq2))
    seq1 = seq1[:min_length]
    seq2 = seq2[:min_length]
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    similarity = matches / min_length
    
    # Apply weighted scoring for consecutive matches
    consecutive_bonus = 0
    current_streak = 0
    
    for a, b in zip(seq1, seq2):
        if a == b:
            current_streak += 1
            consecutive_bonus += (current_streak * 0.1)  # Bonus increases with streak length
        else:
            current_streak = 0
    
    # Normalize bonus and combine with base similarity
    final_similarity = min(1.0, similarity + (consecutive_bonus / len(seq1)))
    return final_similarity

def search_gene(self, sequence: str) -> Dict[str, List[Tuple[str, float]]]:
    """
    Search for matching gene sequences and return matches with similarity scores.
    Added to FusionTree class.
    """
    matches = {}
    
    for disease, gene_sequences in self.genes.items():
        disease_matches = []
        
        for gene in gene_sequences:
            similarity = calculate_sequence_similarity(sequence, gene)
            if similarity > 0.6:  # Threshold for considering a match
                disease_matches.append((gene, similarity))
        
        if disease_matches:
            matches[disease] = disease_matches
    
    return matches

if __name__ == "__main__":
    main()