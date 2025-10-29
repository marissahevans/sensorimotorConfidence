# sensorimotorConfidence
Data and code for 'Feeling lucky? Prospective and retrospective cues for sensorimotor confidence'

### Abstract
On a daily basis, humans interact with the outside world using judgments of sensorimotor confidence, constantly evaluating our actions for success. We ask, what sensory and motor-execution cues are used in making these judgements and when are they available? Two sources of temporally distinct information are prospective cues, available prior to the action (e.g., knowledge of motor noise and past performance), and retrospective cues specific to the action itself (e.g., proprioceptive measurements). We investigated the use of these two cues in two tasks, a secondary motor-awareness task and a main task in which participants reached toward a visual target with an unseen hand and then made a continuous judgment of confidence about the success of the reach. Confidence was reported by setting the size of a circle centered on the reach-target location, where a larger circle reflects lower confidence. Points were awarded if the confidence circle enclosed the true endpoint, with fewer points returned for larger circles. This incentivized accurate reaches and attentive reporting to maximize the score. We compared three Bayesian-inference models of sensorimotor confidence based on either prospective cues, retrospective cues, or both sources of information to maximize expected gain (i.e., an ideal-performance model). Our findings primarily showed two distinct strategies: participants either performed as ideal observers, using both prospective and retrospective cues to make the confidence judgment, or relied solely on prospective information, ignoring retrospective cues. Thus, participants can make use of retrospective cues, evidenced by the behavior observed in our motor-awareness task, but these cues are not always included in the computation of sensorimotor confidence.

### Participants
16 right-handed participants were recruited from the New York University student population (mean age = 25.5 years, SD = 3 years, two male). All but one participant were naïve to the design of the experiment and only naïve participants are presented as examples in this paper. All participants had normal or corrected-to-normal vision, no physical limitations of the right arm and no self-reported motor abnormalities. All participants were tested in both the control motor-awareness task (Task 1) and the main confidence-judgment task (Task 2).

### Experimental Files
The experimental set up involves a Wacom Cintiq 22 tablet placed on a table below a projector. The experimental code runs under the assumption there are 3 screens in use (main computer, tablet and projector). PsychToolBox runs all visualizations through Matlab

Data is collected through the pen position via GetMouse() and WinTabMex(), and with a Griffin PowerMate control knob using PsychPowerMate()

'metaReachExpCode.m' is the main shell for running the experiment. This is where you can enter the participant specifications and request the control experiment or practice trials. This calls all the other functions in the file accordingly.

### Analysis
The primary analysis was performed using the fmincon.m function. Model recovery run with parallel processing on an HPC cluster. All relevant files are in the 'data' folder but to follow the process yourself run the files in this order: 
1: lookUpTabGenModelX.m (this generates the lookup-tables for maximum expected gain based on a set of paramters for each model.) 
2: mmTransform.m (puts the raw data into mm coordinates and centers all trials)
3: participantDataFit.m (runs fmincon.m to minimize negative log liklihood fit)
4. dataSimulatorModelOnly.m (simulates data sets from the models)
5. fminConModelXparfor300.m (performs the model recovery on the simulated data to determin the original model used to create them)
