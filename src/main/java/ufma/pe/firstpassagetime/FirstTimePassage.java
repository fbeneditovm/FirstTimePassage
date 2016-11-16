/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package ufma.pe.firstpassagetime;

import java.math.BigDecimal;
import java.math.RoundingMode;

/**
 *
 * @author fbeneditovm
 */
public class FirstTimePassage {
    
    //The transition Matrix
    BigDecimal[][] tMatrix;
    
    
    /**
     * Constructor
     * @param matrix the transition tMatrix
     */
    public FirstTimePassage(BigDecimal[][] matrix){
        this.tMatrix = matrix;
    }
    
    /**
     * Calculates a specific probability (Fij)^n
     * which is the chance of reachin J from I
     * after n iterations
     * @param i the initial state
     * @param j the desired state
     * @param n the number of iterations
     * @return (Fij)^n
     */
    public BigDecimal probOfReachJfromI(int i, int j, int n){
        
        if(n==1)
            return tMatrix[i][j];  
        
        BigDecimal P = new BigDecimal("0");
        
        for(int k=0; k<tMatrix.length; k++){
            P = P.add(tMatrix[i][k].multiply(probOfReachJfromI(k, j, n-1)));
        }
        return P;
    }
    
    /**
     * Calculates the Stead State Probabilities
     * @return an array with the Stead State Probabilities
     */
    public BigDecimal[] steadStateProb(){
        
        int size = tMatrix.length;
        
        //A is the equations' tMatrix of coeficients
        BigDecimal[][] matrixA = new BigDecimal[size][size];
        
        //The equations' vector of coeficients
        BigDecimal[] vectorB = new BigDecimal[size];
        
        /**
         * The next two loops construct the equations
         */
        for(int l=0; l<(size); l++){
            //The independent term is aways 0 based on the formula
            vectorB[l]= new BigDecimal("0");
            
            for(int k=0; k<(size); k++){
                //atribution of the coeficients based on the formula
                /*
                matrixA[l][k] = (l==k)
                        ? tMatrix[k][l].subtract(new BigDecimal("1.0"))
                        : tMatrix[k][l];
                */
                if(l==k)
                    matrixA[l][k] = tMatrix[k][l].subtract(new BigDecimal("1.0"));
                else
                    matrixA[l][k] = tMatrix[k][l];
            }
        }
        BigDecimal[] incognito = GaussianElimitationBD.lsolve(matrixA, vectorB);
        
        return incognito;
    }
    
     /**
     * Calculates a particular Stead State Probability
     * @return  the Stead State Probability
     */
    public BigDecimal steadStateProb(int state){
        
        int size = tMatrix.length;
        
        //A is the equations' tMatrix of coeficients
        BigDecimal[][] matrixA = new BigDecimal[size][size];
        
        //The equations' vector of coeficients
        BigDecimal[] vectorB = new BigDecimal[size];
        
        /**
         * The next two loops construct the equations
         */
        for(int l=0; l<(size); l++){
            //The independent term is aways 0 based on the formula
            vectorB[l]= new BigDecimal("0");
            
            for(int k=0; k<(size); k++){
                //atribution of the coeficients based on the formula
                /*
                matrixA[l][k] = (l==k)
                        ? tMatrix[k][l].subtract(new BigDecimal("1.0"))
                        : tMatrix[k][l];
                */
                if(l==k)
                    matrixA[l][k] = tMatrix[k][l].subtract(new BigDecimal("1.0"));
                else
                    matrixA[l][k] = tMatrix[k][l];
            }
        }
        BigDecimal[] incognito = GaussianElimitationBD.lsolve(matrixA, vectorB);
        
        return incognito[state];
    }
    
    /**
     * Calculates a Mean First Time Passage Uij
     * @param i the initial state
     * @param j the desired state
     * @return u
     */
    public BigDecimal meanFirstTimePassage(int i, int j){
        
        BigDecimal result = new BigDecimal("1.0");
        if(i==j){
            result = result.divide(steadStateProb(j), 8, RoundingMode.HALF_UP);
            return result;
        }
        
        int size = (tMatrix.length-1);
        //A is the equations' tMatrix of coeficients
        BigDecimal[][] matrixA = new BigDecimal[size][size];
        
        //The equations' vector of coeficients
        BigDecimal[] vectorB = new BigDecimal[size];
        
        /**
         * The next two loops construct the equations
         * based on j and the Transition Matrix
         */
        for(int x=0; x<(size); x++){
            //The independent term is aways -1 based on the formula
            vectorB[x]= new BigDecimal("-1");
            //x is the actual equation and l is the actual line in tMatrix
            //the j line is skipped
            int l = (x<j) ? x : (x+1);
            for(int y=0; y<(size); y++){
                //y is the actual coeficient of the equation and k is the column in tMatrix
                //the j column is skipped
                int k = (y<j) ? y : (y+1);
                //atribution of the coeficients based on the formula
                /*
                matrixA[x][y] = (l==k)
                        ? tMatrix[l][k].subtract(new BigDecimal("1.0"))
                        : tMatrix[l][k];
                */
                if(l==k)
                    matrixA[x][y] = tMatrix[l][k].subtract(new BigDecimal("1.0"));
                else
                    matrixA[x][y] = tMatrix[l][k];
            }
        }
        BigDecimal[] incognito = GaussianElimitationBD.lsolve(matrixA, vectorB);
        /*
         * incognito contains the values of all Ukj.
         * Since the k skips j,
         * the desired Mean First Time Passage Uij
         * is the incognito[i] if i<j or
         * incognito[i-1] if i>j
         */
        result = (i<j) ? incognito[i] : incognito[i-1];
        
        return result;
    }
}
