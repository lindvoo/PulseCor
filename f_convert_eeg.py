
import mne

# Convert EEG files
def readineeg(filename,dSF):

  raw = mne.io.read_raw_brainvision(filename)
  raw.load_data()  # preload so you can resample
  raw.resample(dSF)  # resample to 50 HZ
  dat = raw.get_data()  # extract data

  # Save each channel in a txt file
  for c, val in enumerate(raw.ch_names):
      thefile = open(filename[:-5] + "_" + val + ".txt", 'w')
      for item in dat[c]:
          thefile.write("%s\n" % item)
      thefile.close()

  # Get stimulus events
  events = mne.events_from_annotations(raw)
  thefile = open(filename[0][:-5] + "_stimcodes.txt", 'w')
  for c, val in enumerate(list(events[0])):
      if val[0] > 0:
          thefile.write("%s\n" % val[0])
  thefile.close()
  print("files coverted!")
