# Costas Array Generator

This repository contains programs that generate Costas arrays using the following fundamental algebraic methods:
1. Basic Methods:
   
    •	Welch
  
    •	Lempel
  
    •	Golomb
2. Variations:
   
    •	Welch-0
  
    •	Taylor-4
  
    •	Taylor-T1-T0
  
    •	Golomb-G4-G5
  
    •	Polymorphs:
      Welch, Lempel, Golomb
   
## **Theoretical Background**

For detailed information on the methods, you can refer to my paper and MSc thesis:

•	F. Güleç and E. Afacan, "Usage of Costas arrays in Low Probability of Intercept radars," 2015 23nd Signal Processing and Communications Applications Conference (SIU), Malatya, Turkey, 2015, pp. 347-350, https://ieeexplore.ieee.org/abstract/document/7129830.

•	F. Gulec, "Usage of Costas arrays in Low Probability of Intercept radars," MSc Thesis, Gazi University, 2015. Link: https://tez.yok.gov.tr/UlusalTezMerkezi/tezDetay.jsp?id=BV5vRnsGI968q47wRc2_HQ&no=Ej9q9ft9sPwGMpdRq8uB8g

While these documents are in Turkish, the methods used are well-documented in the literature, such as in works by Drakakis or Golomb.

## **Key Details**

•	These methods generate Costas arrays for prime or prime powers.
•	Variations and polymorphs may only support a limited range of orders.
•	Welch is the fastest method, followed by Lempel. Golomb and its variations may take longer due to their complexity.

## **Usage**

Outputs are saved as .txt files in the working directory. You can input the order of the Costas array when prompted in the terminal. Outputs are written to .txt files in the same folder. After inputting the desired order of Costas arrays in the terminal, the program will print the Costas arrays and the primitive elements. You will find some example results in the folders.
The code is written in C++ (mostly C-style syntax). For Lempel and Golomb-based methods, a parameter N should be adjusted for orders greater than 100 to prevent anomalous results.

I wrote these codes for my MSc thesis in 2014, and I am sharing them 10 years later. Contributions are welcome!

## **Contributions**

 If you use or improve upon this code, please cite my work as follows:



## **Citation Requirement**

If you use or build upon this code in your research, or if this code is used for any academic work, publication, or research, proper attribution and citation of the paper is **required**. Please cite the paper in any related publications, presentations, or derived works.

•	F. Güleç and E. Afacan, "Usage of Costas arrays in Low Probability of Intercept radars," 2015 23nd Signal Processing and Communications Applications Conference (SIU), Malatya, Turkey, 2015, pp. 347-350, doi: 10.1109/SIU.2015.7129830. https://ieeexplore.ieee.org/abstract/document/7129830
