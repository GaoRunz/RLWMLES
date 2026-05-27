If you would like to run RL-WMLES, you first need to install the Smarties library. Detailed installation instructions can be found at:
https://github.com/cselab/smarties

Next, you need to couple our WMLES code with Smarties. This can be done by simply adding the provided code to the `apps` directory of Smarties. In this framework, our WMLES solver serves as the reinforcement learning environment interacting with the RL agent. The remaining training and testing procedures are consistent with those described in the Smarties documentation.

Finally, modifications related to the WMLES solver can be found in:
`RLcode/apps/WMLES_1`
