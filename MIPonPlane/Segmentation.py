import math
import numpy as np
import SimpleITK as sitk
import matplotlib.pyplot as plt

from matplotlib.path import Path
plt.ion()


class GMM_Segmentation (object):
    """
    img = np array, nCluster = n classi
    seg_img?????? output

    """
    def __init__(self, img, nCluster, likelihoodTreshold=None, maxIter= None, beta = None, nNeighbors = None):
        self.img = img #np array
        self.nCluster = nCluster #int, numb classes
        self.beta = beta #MRF parameter
        self.nNeighbors = nNeighbors #4 or 8, MRF paraneter
        self.MRFflag = False #no MRF computation
        self.exit_conditions = False

        #likelihood and n max iteration initialization in case of None
        if likelihoodTreshold is None:
            self.likelihoodTreshold = 130
            print('Absent loglikelihood value, computation with default value'), self.likelihoodTreshold
        else:
            self.likelihoodTreshold = likelihoodTreshold
        if maxIter is None:
            self.maxIter = 300
            print('Absent number max iterations, computation with default value'), self.maxIter
        else:
            self.maxIter = maxIter
        if (self.beta is not None or self.nNeighbors is not None):
            self.MRFflag = True #avvio MRF in caso di chiamata con almeno uno dei 2 parametri caratteristici
            if self.beta is None:
                self.beta = 1.5
                print('Absent beta value, computation with default value'), self.beta
            else:
                self.beta = beta
            if self.nNeighbors is None:
                self.nNeighbors = 4
                print('Absent number of neighbors, computation with default value'), self.nNeighbors
            else:
                self.nNeighbors = nNeighbors


        #initialization to enter in the first Exp-Max step
        self.loglikelihood = 1000
        self.diff_loglikelihood = 20000
        self.iter = 0
        self.seg_img = img
        self.inMRF = False

        #IMAGES PARAMETERS
        self.minInt = np.amin(self.img)
        self.maxInt = np.amax(self.img)
        self.nRow = self.img.shape[0]
        self.nCol = self.img.shape[1]
        self.nPixel = self.nRow * self.nCol

        #INITIALIZATION (mean, variance)
        #mu = (1 : nCluster) * maxVal / (nCluster + 1);
        self.mu = (np.arange(self.nCluster) + np.ones(self.nCluster)) * (self.maxInt / (self.nCluster + 1))
        #vari = ones(1, nCluster) * maxVal;
        self.vari = np.ones([self.nCluster]) * self.maxInt

        # INITIALIZATION (probmatrix, sum of each pixel distr ~ sum_gauss, proportion of each cluster ~ pr)
        self.prob_matrix = np.zeros([self.nRow, self.nRow, self.nCluster])
        self.sum_gauss = np.zeros([self.nRow, self.nCol])
        self.pr = np.ones([self.nCluster]) / self.nCluster
        self.distr = np.zeros([self.nRow, self.nCol, self.nCluster])

    #CALLING HERE THE PROPER FUNCTION
    def compute(self):
        self.seg_img = self.expectation_maximization_steps()
        if self.MRFflag is True:
            self.inMRF = True
            self.exit_conditions = False
            self.loglikelihood = 20000
            self.seg_img = self.expectation_maximization_steps()
            return self.seg_img
        else:
            return self.seg_img

    #EXPECTATION MAXIMIZATION STEPS (BOTH SIMPLE ~inMRF=False~ AND MRF ~inMRF=True)
    def expectation_maximization_steps (self):
        while(self.exit_conditions == False):
            self.iter += 1
            if self.inMRF is False:
                print self.iter, '. Iteration (EM)'
            else:
                print self.iter, '. Iteration (EM - MRF)'
                self.matrix_prior = self.matrix_prior_computation()

            #EXPECTATION STEP
            loglikelihood_old = self.loglikelihood
            self.distr = self.distr_computation()
            self.sum_gauss = self.sum_gauss_computation()
            self.loglikelihood = self.loglikelihood_computation()
            self.prob_matrix = self.prob_matrix_computation()
            self.diff_loglikelihood = abs(loglikelihood_old - self.loglikelihood)
            print 'Diff_loglikelihood', self.diff_loglikelihood
            #MAXIMIZATION STEP
            self.pr = self.pr_computation()
            self.mu = self.mu_comp()
            self.vari = self.vari_comp()
            #UPDATE EXIT CONDITIONS
            self.exit_conditions = self.exit_conditions_computation()
        else:
            if self.iter == self.maxIter:
                print('ERROR: Max iter')
                return -1
            return self.seg_img

    def exit_conditions_computation(self):   #VERIFY EXIT CONDITIONS
        if self.diff_loglikelihood <= self.likelihoodTreshold:
            print '\nConvergenza raggiunta in ', self.iter, ' iterazioni\n'
            self.seg_img = np.argmax(self.prob_matrix, 2)
            self.exit_conditions = True
        elif self.iter == self.maxIter:
            #display di errore per maxIter raggiunto
            self.exit_conditions = True
        else:
            self.exit_conditions = False
        return self.exit_conditions

    def distr_computation(self):
        #costn(l, c, z) = pr(z) / sqrt(2 * pi * vari(z));
        #espon(l, c, z) = -((img(l, c) - mu(z)). ^ 2 / (2 * vari(z)));
        #distr(l, c, z) = costn(l, c, z) * exp(espon(l, c, z));
        costn = np.zeros([self.nRow, self.nCol, self.nCluster])
        espon = np.zeros([self.nRow, self.nCol, self.nCluster])

        for l in range(self.nRow):
            for c in range(self.nCol):
                for z in range(self.nCluster):
                    if self.inMRF is False:
                        costn[l, c, z] = self.pr[z] / math.sqrt(2 * math.pi * self.vari[z])
                    else:
                        costn[l, c, z] = self.matrix_prior[l, c, z] / math.sqrt(2 * math.pi *self.vari[z])
                    espon[l, c, z] = -((self.img[l, c] - self.mu[z])**2 / (2 * self.vari[z]))
                    self.distr[l, c, z] = costn[l, c, z] * np.exp(espon[l, c, z])
        return self.distr

    def sum_gauss_computation(self):
        # somma lungo la terza direzione (0 righe, 1 colonne, 2 terza dim)
        self.sum_gauss = np.sum(self.distr, axis = 2)
        #sum_gauss: nR * nC
        return self.sum_gauss

    def prob_matrix_computation(self):
        prob_matrix = np.zeros([self.nRow, self.nCol, self.nCluster])
        #self.sum_gauss = self.sum_gauss_computation()
        for l in range(self.nRow):
            for c in range(self.nCol):
                for z in range(self.nCluster):
                    #prob_matrix(l, c, z) = distr(l, c, z). / sum_gauss(l, c);
                    prob_matrix[l, c, z] = self.distr[l, c, z] / self.sum_gauss[l, c]
        return prob_matrix

    def loglikelihood_computation(self):
        return (np.sum(np.sum(np.log(self.sum_gauss))))

    def pr_computation(self):
        pr_sum_lines = np.sum(self.prob_matrix_computation(), 0)
        pr_sum = np.sum(pr_sum_lines, 0)
        pr = pr_sum / self.nPixel
        pr = pr / np.sum(pr)
        return pr

    def mu_comp(self):
        s_row = np.sum(self.prob_matrix, 0)
        s_row_line = np.sum(s_row, 0)
        for cl in range(self.nCluster):
            self.mu[cl] = np.sum(np.sum(self.img[:, :] * self.prob_matrix[:, :, cl])) / s_row_line[cl]
        return self.mu

    def vari_comp(self):
        s_row = np.sum(self.prob_matrix, 0)
        s_row_line = np.sum(s_row, 0)
        for cl in range(self.nCluster):
            stnde = np.zeros(self.img.shape)
            stnde[:, :] = self.img[:, :] - self.mu[cl]
            self.vari[cl] = np.sum(np.sum(stnde[:, :]**2 * self.prob_matrix[:, :, cl])) / s_row_line[cl]
        return self.vari

    #da chiamare all'inizio o a exp-max completata
    def original_intensities(self):
        x_intensities = np.arange(self.minInt, self.maxInt + 1)
        y_intensities = np.zeros(x_intensities.shape)
        for intensity in range(self.minInt, self.maxInt + 1):
            if intensity != 0:
                array = (self.img == intensity)
                y_intensities[intensity - self.minInt] = sum(sum(array))
        #x_intensities contiene i valori di ogni intensita' nel range dell'immagine
        #y_intensities[i] contiene il numero di pixel per ogni intensita'
        y_intensities = y_intensities / self.nPixel
        #y_intensities contiene ora la frequenza di ogni intensita'
        return [x_intensities, y_intensities]

    #MRF only
    def matrix_prior_computation(self):
        self.matrix_prior = np.zeros([self.nRow, self.nCol, self.nCluster])
        for row in range(self.nRow):
            for col in range(self.nCol):
                neighborhood = self.find_neighborhood(row, col)
                cl_frequence = np.zeros(self.nCluster)
                for nb in range(neighborhood.size):
                    for cl in range(self.nCluster):
                        if neighborhood[nb] == cl:
                            cl_frequence[cl] += 1
                energy_fun = self.beta * cl_frequence
                self.matrix_prior[row, col, :] = np.exp(energy_fun) / np.sum(np.exp(energy_fun))
        return self.matrix_prior

    #MRF only
    def find_neighborhood(self, row, col):
        upleft = (row == 0 and col == 0)
        upright = (row == 0 and col == self.nCol - 1)
        downright = (row == self.nRow - 1 and col == 0)
        downleft = (row == self.nRow - 1 and col == self.nCol - 1)
        upline = (row == 0)
        downline = (row == self.nRow - 1)
        leftline = (col == 0)
        rightline = (col == self.nCol - 1)
        if (self.nNeighbors == 4 or self.nNeighbors ==8):
            if (upleft or upright or downright or downleft):
                neighborhood = np.zeros(2)
                if upleft:
                    neighborhood[0] = self.seg_img[0, 1]
                    neighborhood[1] = self.seg_img[1, 0]
                elif upright:
                    neighborhood[0] = self.seg_img[0, self.nCol - 2]
                    neighborhood[1] = self.seg_img[1, self.nCol - 1]
                elif downleft:
                    neighborhood[0] = self.seg_img[self.nRow - 2, 0]
                    neighborhood[1] = self.seg_img[self.nRow - 1, 1]
                elif downright:
                    neighborhood[0] = self.seg_img[self.nRow - 2, self.nCol - 1]
                    neighborhood[1] = self.seg_img[self.nRow - 1, self.nCol - 2]
            elif (upline, downline, leftline, rightline):
                neighborhood = np.zeros(3)
                if upline:
                    neighborhood[0] = self.seg_img[0, col - 1]
                    neighborhood[1] = self.seg_img[0, col + 1]
                    neighborhood[2] = self.seg_img[1, col]
                elif downline:
                    neighborhood[0] = self.seg_img[self.nRow - 1, col - 1]
                    neighborhood[1] = self.seg_img[self.nRow - 1, col + 1]
                    neighborhood[2] = self.seg_img[self.nRow - 2, col]
                elif leftline:
                    neighborhood[0] = self.seg_img[row + 1, 0]
                    neighborhood[1] = self.seg_img[row - 1, 0]
                    neighborhood[2] = self.seg_img[row, 1]
                elif rightline:
                    neighborhood[0] = self.seg_img[row + 1, self.nCol - 1]
                    neighborhood[1] = self.seg_img[row - 1, self.nCol - 1]
                    neighborhood[2] = self.seg_img[row, self.nCol - 2]
            else:
                neighborhood = np.zeros(4)
                neighborhood[0] = self.seg_img(row, col + 1)
                neighborhood[1] = self.seg_img(row, col - 1)
                neighborhood[2] = self.seg_img(row + 1, col)
                neighborhood[3] = self.seg_img(row - 1, col)

            if self.nNeighbors == 8:
                if upleft:
                    extra = np.array(self.seg_img[1, 1])
                elif upright:
                    extra = np.array(self.seg_img[1, self.nCol - 2])
                elif downleft:
                    extra = np.array(self.seg_img[self.nRow - 2, 1])
                elif downright:
                    extra = np.array(self.seg_img[self.nRow - 2, self.nCol - 2])
                elif upline:
                    extra = np.array([self.seg_img[1, col + 1], self.seg_img[1, col - 1]])
                elif downline:
                    extra = np.array([self.seg_img[self.nRow - 2, col + 1], self.seg_img[self.nRow - 2, col - 1]])
                elif leftline:
                    extra = np.array([self.seg_img[row + 1, 1], self.seg_img[row + 1, - 1]])
                elif rightline:
                    extra = np.array([self.seg_img[row + 1, self.nCol - 2], self.seg_img[row - 1, self.nCol - 2]])
                else:
                    extra = np.array([self.seg_img[row + 1, col + 1], self.seg_img[row + 1, col - 1], self.seg_img[row - 1, col + 1], self.seg_img[row - 1, col -1]])
                neighborhood = np.append(extra, neighborhood)

        return neighborhood