# LEO Satellite Diversity
LEO Satellite Diversity is a small project for researching the effects of spatial diversity when
applied on a non-terrestrial network, and specifically
a satellite network.

## Prerequisites
- MATLAB R2024b or later
- MATLAB Satellite Toolbox
- MATLAB Antenna Toolbox
- MATLAB Phased Array Toolbox

## Simulation Scenarios
### Simple MIMO 2x2
**System Model**: Two Ground Stations (GS) and two Satellites (SAT). Each ground station has two gimbals with one transceiver each. Each satellite gimbal points at one satellite. Thus, two uncorrelated paths are established:
- Path 1: *GS1Tx1 == Sat1Rx == Sat1Tx == GS2Rx1*
- Path 2: *GS1Tx2 == Sat2Rx == Sat2Tx == GS2Rx2*

The received ${E_b/N_0}$ is calculated for the whole multi-hop link. The model is evaluated using the BER metric for different digital modulations, i.e. QPSK, 8-QAM, 16-QAM.

### MegaLEO Constellation
**System Model**: Two GS and 1000 satellites. During each communication occassion the transmit GS selects the satellite closer to the target GS. The satellite does the same, and the routing is performed. When the route is established, the device chain is used to create a cell array of link objects. The BER of the link chain is calculated.

### Adaptive Combining Scheme
**Motivation**: Channels are not always known to the receiver. The receiver estimates the channel by using pilot bits (pilot SNR). The estimation can be lossy, meaning that MRC (maximal ratio combining) can produce low quality results. SC (selection combining) is simpler, and requires less resources than MRC, if the channel estimation is uncertain, better use SC instead of MRC.

**System Model**: Simple MIMO 2x2 is used as the main system model. Assuming optimal alignment of the received symbols in the combiner, MRC produces the optimal result but with a possibly uncertain channel estimate.

**Uncertainty Model**: Modelled as a ratio of the estimated channel coefficient $\hat{h}$ and the estimation error variance ${\sigma_e^2}$. If the ratio is close to 1, then the estimation is equal to the estimation noise, and the coefficient is useless. If the ratio is close to 0, then the estimation can be used. A threshold is defined as ${\gamma}$, where if it is 1, then only MRC is used, and if it is 0, only SC is used. The switching occurs based on ${\gamma}$.

### Adaptive Combining Scheme with Reinforcement Learning

**Motivation**: Using a hard-coded threshold is pointless, thus a Q-Learning Agent is trained to select the combining technique. The agent creates a Q-Matrix, approximating the true Quality Function of each action. The discrete nature of the Q-Matrix limits the precision of the reward. After the Q-Matrix creation, the model is deployed on random SNR values between -10 and 20 dB. The results display that the agent is capable of distinguishing between low and high uncertainty and selecting the best combining technique each time.
