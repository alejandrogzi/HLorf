# TRANSAID: A Hybrid Deep Learning Framework for Translation Site Prediction

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**TRANSAID (TRANSlation AI for Detection)** is a powerful and accurate tool for identifying translation initiation sites (TIS) and translation termination sites (TTS) from full-length RNA transcripts. It leverages a novel hybrid deep learning architecture, combining dilated convolutions and residual blocks, to predict complete Open Reading Frames (ORFs) with high precision.

This tool was designed to overcome the limitations of existing methods by being trained on both protein-coding (NM) and non-coding (NR) transcripts, significantly reducing false-positive predictions on non-coding RNAs while maintaining high sensitivity for genuine ORFs.

## Key Features

-   **High Accuracy**: State-of-the-art performance in predicting perfect ORFs on coding transcripts and correctly identifying non-coding transcripts.
-   **End-to-End Prediction**: Directly predicts complete TIS-TTS pairs from raw transcript sequences.
-   **Integrated Biological Scoring**: Incorporates a sophisticated scoring system including Kozak context, Codon Adaptation Index (CAI), and GC content to enhance prediction reliability.
-   **Translation Product Generation**: Automatically translates predicted ORFs into protein sequences (`.faa` file).
-   **Cross-Species Applicability**: While trained on human data, the model demonstrates strong performance across diverse eukaryotic species.
-   **GPU Accelerated**: Utilizes GPU for fast predictions on large datasets.

## Online Web Server

For quick and easy analysis of a few sequences without any installation, we provide a user-friendly web server.

-   **URL**: [http://58.242.248.157:6005](http://58.242.248.157:6005)  
-   **Functionality**:
    -   Paste FASTA sequences directly or upload a FASTA file.
    -   Adjust key prediction parameters (`filter_mode`, `integrated_cutoff`, etc.).
    -   Download results as CSV and protein sequences as FAA files.

## Installation

TRANSAID is a command-line tool written in Python and requires a Conda environment with specific dependencies, including PyTorch with GPU support.

### 1. Prerequisites

-   **Anaconda or Miniconda**: To manage the environment.
-   **NVIDIA GPU**: With CUDA drivers installed (version 11.x or 12.x recommended).

### 2. Installation Steps

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/your-username/transaid.git
    cd transaid
    ```

2.  **Create and activate the Conda environment:**
    We provide a `requirements.txt` file. You can create a new environment and install the packages.

    ```bash
    # Create a new conda environment with Python 3.11
    conda create -n transaid_env python=3.11 -y

    # Activate the environment
    conda activate transaid_env

    # Install PyTorch with CUDA support (IMPORTANT: choose the command that matches your CUDA version)
    # For CUDA 12.1:
    pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121
    # For CUDA 11.8:
    # pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

    # Install other required packages
    pip install -r requirements.txt
    ```
    Your `requirements.txt` should contain at least:
    ```
    numpy
    pandas
    biopython
    tqdm
    ```

3.  **Download Pre-trained Models:**
    Download the pre-trained model files and place them in the `model/` directory.
    -   `TRANSAID_Embedding_batch4_best_model.pth`: Trained on a mixed dataset of NM and NR transcripts (Recommended).
    -   `TRANSAID_Embedding_batch4_TrainOnlyNM_best_model.pth`: Trained on NM transcripts only (for comparison).

## Command-Line Usage

The main script for making predictions is `predict.py`.

### Basic Example

Here is a basic command to predict ORFs on a FASTA file:

```bash
python predict.py \
    --input /path/to/your/transcripts.fna \
    --model_path ./model/TRANSAID_Embedding_batch4_best_model.pth \
    --output /path/to/your/output_prefix \
    --gpu 0
```

This will generate two files:
-   `/path/to/your/output_prefix.csv`: A detailed CSV file with all predicted ORFs and their scores.
-   `/path/to/your/output_prefix.faa`: A FASTA file containing the translated protein sequences of the predicted ORFs that passed the filter.

### Parameter Explanations

Here is a detailed explanation of all available command-line options:

| Argument | Description | Default |
| :--- | :--- | :--- |
| `--input` | **[Required]** Path to the input FASTA file containing RNA transcript sequences. | |
| `--model_path`| **[Required]** Path to the pre-trained `.pth` model file. | |
| `--output` | **[Required]** The base path and prefix for the output files (e.g., `results/my_preds`). The script will append `.csv` and `.faa`.| |
| `--gpu` | GPU device ID to use. Set to `-1` to run on CPU (slower). | `0` |
| `--filter_mode` | Prediction filtering mode. `best`: keeps only the highest-scoring ORF per transcript. `all`: keeps all ORFs that pass the score cutoffs. | `best` |
| `--integrated_cutoff` | The minimum `Integrated_Score` an ORF must have to be considered a positive prediction. This is the primary filtering threshold. | `0.635` |
| `--orf_length_cutoff`| The minimum length (in amino acids) of the translated protein for an ORF to be reported. | `30` |
| `--batch_size`| The number of sequences to process in one batch. Adjust based on GPU memory. | `4` |
| `--max_seq_len` | The maximum sequence length the model was trained on. Sequences longer than this will be truncated. | `27112` |
| `--kozak_cutoff`| Minimum Kozak sequence score for filtering. Only applied if `integrated_cutoff` is not the main focus. | `0.05` |
| `--tis_cutoff`| Minimum TIS probability score from the model. | `0.1` |
| `--tts_cutoff`| Minimum TTS probability score from the model. | `0.1` |
| `--save_raw_predictions`| If specified, saves the raw model output probabilities as a `.pkl` file for debugging. | `False` |

---

## Citing TRANSAID

If you use TRANSAID in your research, please cite our paper:

> [Wu, Z., Wang, B., et al. (2025). TRANSAID: A Hybrid Deep Learning Framework... ]

We hope TRANSAID accelerates your research in transcriptomics and proteomics!
