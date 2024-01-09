class Dna {
  static ACGT = new Map([
    ["00", "A"],
    ["01", "C"],
    ["10", "G"],
    ["11", "T"],
  ]);

  static AMINO_ACIDES_CODONS = {
    A: ["GCU", "GCC", "GCA", "GCG"],
    R: ["CGU", "CGC", "CGA", "CGG"],
    N: ["AAU", "AAC"],
    D: ["GAU", "GAC"],
    C: ["UGU", "UGC"],
    Q: ["CAA", "CAG"],
    E: ["GAA", "GAG"],
    G: ["GGU", "GGC", "GGA", "GGG"],
    H: ["CAU", "CAC"],
    I: ["AUU", "AUC", "AUA"],
    L: ["CUU", "CUC", "CUA", "CUG"],
    K: ["AAA", "AAG"],
    M: ["AUG"],
    F: ["UUU", "UUC"],
    P: ["CCU", "CCC", "CCA", "CCG"],
    S: ["UCU", "UCC", "UCA", "UCG"],
    T: ["ACU", "ACC", "ACA", "ACG"],
    W: ["UGG"],
    Y: ["UAU"],
    V: ["GUU", "GUC", "GUA", "GUG"],
    B: ["UAA", "UGA", "UAG"],
    O: ["UUA", "UUG"],
    U: ["AGA", "AGG"],
    X: ["AGU", "AGC"],
    Z: ["UAC"],
  };

  static AminoAcidsToDNA(plainText, ambiguity) {
    let dnaSequence = "";

    for (let i = 0; i < plainText.length; i++) {
      const aminoAcid = plainText[i];
      const currentAmbiguity = ambiguity[i];

      const dnaChar = this.AMINO_ACIDES_CODONS[aminoAcid][currentAmbiguity];

      dnaSequence += dnaChar;
    }

    return dnaSequence;
  }

  static DNAToAminoAcids(dna) {
    let aaString = "";
    let ambiguityString = "";

    const fillerCount = dna.length % 3 > 0 ? 3 - (dna.length % 3) : 0;

    dna += Array(fillerCount).fill("G").join("");

    for (let i = 0; i < dna.length; i += 3) {
      const currentCodon = dna.slice(i, i + 3);
      const foundCodon = Object.keys(this.AMINO_ACIDES_CODONS).find((key) =>
        this.AMINO_ACIDES_CODONS[key].includes(currentCodon)
      );

      if (foundCodon) {
        const index =
          this.AMINO_ACIDES_CODONS[foundCodon].indexOf(currentCodon);
        ambiguityString += index;

        aaString += foundCodon;
      }
    }

    return [aaString, ambiguityString];
  }

  static DNAToRNA(dna) {
    return dna.replaceAll("U", "T");
  }

  static DNAToBinary(dna) {
    const rna = this.DNAToRNA(dna);
    let binary = "";

    for (let i = 0; i < rna.length; i++) {
      const bin = [...this.ACGT].find(([_, value]) => rna[i] === value)[0];
      binary += bin;
    }
    return binary;
  }

  static RNAtoDNA(rna) {
    return rna.replaceAll("T", "U");
  }

  static binaryToDNA(bin) {
    let dnaString = "";

    for (let i = 0; i < bin.length; i += 2) {
      const adn = this.ACGT.get(bin[i] + bin[i + 1]);
      dnaString += adn;
    }

    dnaString = this.RNAtoDNA(dnaString);

    return dnaString;
  }
}

class MathManipulator {
  static modInverse(a, m) {
    for (let i = 0; i < m; i++) {
      if ((a * i) % m === 1) {
        return i;
      }
    }
    return 1;
  }

  static determinant(matrix) {
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
  }

  static inverseMatrix(matrix) {
    const det = this.determinant(matrix);
    const modInv = this.modInverse(det, 26);

    const inverseMatrix = [
      [matrix[1][1], -matrix[0][1]],
      [-matrix[1][0], matrix[0][0]],
    ];

    // Multiply each element by the modular inverse of the determinant
    for (let i = 0; i < 2; i++) {
      for (let j = 0; j < 2; j++) {
        inverseMatrix[i][j] = (inverseMatrix[i][j] * modInv) % 26;
        // Ensure the result is non-negative
        if (inverseMatrix[i][j] < 0) {
          inverseMatrix[i][j] += 26;
        }
      }
    }

    return inverseMatrix;
  }
}

class StringManipulator {
  static charToNumber(char) {
    return char.charCodeAt(0) - "A".charCodeAt(0); // 65
  }

  static numberToChar(number) {
    return String.fromCharCode(number + "A".charCodeAt(0)); // 65
  }

  static binaryToText(bin) {
    let result = [];
    for (let i = 0; i < bin.length; i += 8) {
      result.push(bin.slice(i, i + 8));
    }

    return result
      .map((value) => String.fromCharCode(parseInt(value, 2)))
      .join("");
  }

  static stringToBinary(str) {
    let binaryString = "";

    for (let i = 0; i < str.length; i++) {
      const binary = str[i].charCodeAt(0).toString(2).padStart(8, "0");
      binaryString += binary;
    }

    return binaryString;
  }
}

class HillEncrypt {
  static encrypt(plainText, key) {
    let encryptedText = "";

    // append J if odd, we chose J because it doesn't have an amino acidic representation
    if (plainText.length % 2 !== 0) {
      plainText += "J";
    }

    // turn letter into numbers
    for (let i = 0; i < plainText.length; i += 2) {
      // turn letters into numbers
      const letter1 = StringManipulator.charToNumber(plainText[i]);
      const letter2 = StringManipulator.charToNumber(plainText[i + 1]);

      // multiply the lists
      const result1 = (letter1 * key[0][0] + letter2 * key[0][1]) % 26;
      const result2 = (letter1 * key[1][0] + letter2 * key[1][1]) % 26;

      // turn numbers back to letters
      encryptedText +=
        StringManipulator.numberToChar(result1) +
        StringManipulator.numberToChar(result2);
    }

    return encryptedText;
  }
  static textEncrypt(text, KEY) {
    const bin = StringManipulator.stringToBinary(text);
    const dna = Dna.binaryToDNA(bin);
    const [codons, ambiguity] = Dna.DNAToAminoAcids(dna);
    const cipher = this.encrypt(codons, KEY);

    console.log("- Hill DNA Encryption:\n");
    console.log("> Plain text:   ", text, "\n");
    console.log("> binary text:  ", bin, "\n");
    console.log("> DNA sequence: ", dna, "\n");
    console.log("> AA codons:    ", codons, "\n");
    console.log("> Ambiguity:    ", ambiguity, "\n");
    console.log("> Cipher text:  ", cipher, "\n");

    return { cipher, ambiguity };
  }
}

class HillDecrypt {
  static decrypt(cipher, KEY) {
    const keyMatrixInverse = MathManipulator.inverseMatrix(KEY);

    let plainText = "";

    for (let i = 0; i < cipher.length; i += 2) {
      // Convert letters to numbers
      const letter1 = StringManipulator.charToNumber(cipher[i]);
      const letter2 = StringManipulator.charToNumber(cipher[i + 1]);

      // Multiply the inverse key matrix with the column vector of letters
      const result1 =
        (keyMatrixInverse[0][0] * letter1 + keyMatrixInverse[0][1] * letter2) %
        26;
      const result2 =
        (keyMatrixInverse[1][0] * letter1 + keyMatrixInverse[1][1] * letter2) %
        26;

      // Convert the results back to letters
      const decryptedLetter1 = StringManipulator.numberToChar(result1);
      const decryptedLetter2 = StringManipulator.numberToChar(result2);

      // Append the decrypted letters to the plaintext
      plainText += decryptedLetter1 + decryptedLetter2;
    }

    if (plainText.charAt(plainText.length - 1) === "J")
      plainText = plainText.slice(0, -1);

    return plainText;
  }
  static textDecrypt(cipher, ambiguity, KEY) {
    const codons = this.decrypt(cipher, KEY);
    const dnaSequence = Dna.AminoAcidsToDNA(codons, ambiguity);
    const binaryText = Dna.DNAToBinary(dnaSequence);
    const plainText = StringManipulator.binaryToText(binaryText);

    console.log("- Hill DNA Decryption:\n");
    console.log("> Cipher text:  ", cipher, "\n");
    console.log("> Ambiguity:    ", ambiguity, "\n");
    console.log("> AA codons:    ", codons, "\n");
    console.log("> DNA sequence: ", dnaSequence, "\n");
    console.log("> Binary text:  ", binaryText, "\n");
    console.log("> Plain text:   ", plainText, "\n");
  }
}

const KEY = [
  [6, 7],
  [3, 9],
];

const plainText = "Papier extra blanc, qualite superieure";

const { cipher, ambiguity } = HillEncrypt.textEncrypt(plainText, KEY);

console.log("---------------------------------------------------------");

HillDecrypt.textDecrypt(cipher, ambiguity, KEY);
