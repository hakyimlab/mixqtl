# Update note

* v0.2 (8/4/20): Updated the step 1 of multi-snp approaches (mixfine and mixpred). In details, we replace the cv.glmnet based approach to estimate variance by random/mixed effect approach (via emma [source](http://mouse.cs.ucla.edu/emma/index.html)). And we also changed the way of accounting for intercept in TRC part by an easier approach, simply regressing out the intercept from both y and X.
