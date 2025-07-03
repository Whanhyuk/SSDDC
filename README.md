# DNA-Based Data Encoding and Decoding System

This project implements an encoding and decoding system for digital data storage in DNA sequences, featuring robust error correction for insertion, deletion, and substitution errors. The system is based on custom generator matrices and syndrome decoding.

---

## Overview

This Python program enables the transformation of text files into DNA sequences and the accurate retrieval of the original text from DNA-encoded files. It uses error-correcting codes, matrix transformations, and DNA-specific patterns to ensure reliable data storage and recovery.

---

## Key Features

- **Text ↔ DNA Conversion**  
  Converts text to binary, binary to DNA bases, and vice versa, using predefined dictionaries.

- **Matrix-Based Encoding**  
  Utilizes generator matrices to encode data blocks, enabling linear block codes and error correction.

- **Error Correction**  
  Handles both substitution and indel (insertion/deletion) errors with syndrome decoding and coset leaders.

- **Command-Line Interface**  
  Simple command-line arguments for encoding (`-e`) and decoding (`-d`).

---

## How It Works

### Encoding

1. Reads text file input.
2. Converts text to binary, then maps binary to DNA bases (A, T, C, G).
3. Represents DNA bases as matrices for further encoding.
4. Applies generator matrix multiplication to produce codewords.
5. Inserts specific patterns into DNA to help correct deletion errors.
6. Writes the DNA-encoded output to a file.

### Decoding

1. Reads DNA-encoded input file.
2. Detects and corrects deletion errors.
3. Uses syndrome decoding with coset leaders to fix substitution errors.
4. Converts corrected DNA back to binary and then to the original text.
5. Writes the recovered text to a file.

---

## Typical Workflow

### Encoding

```bash
python final_version_pub.py -e input.txt
```

This command will create a file named `DNA_Encoded_input.txt`.

### Decoding

```bash
python final_version_pub.py -d DNA_Encoded_input.txt
```

This will recover the original text as `Retrieved_DNA_Encoded_input.txt`.

---

## Strengths

- **Robust Error Correction**  
  Effective against both deletion and substitution errors, which are crucial for practical DNA storage.

- **Modular and Extensible**  
  Code structure allows for modifications, such as using different generator matrices or error correction strategies.

- **Educational**  
  Clearly demonstrates text-to-binary, binary-to-DNA, and matrix-based error correction processes.

---

## Applications

- DNA data storage research and prototyping
- Teaching concepts of coding theory, error correction, and bioinformatics
- Demonstrating digital-to-biological information encoding

---

## License

This project is released under the MIT License.

---

## Contact

For questions or contributions, please open an issue or submit a pull request.
