
import bioread

# convert ACQ files
def readinacq(filename,dSF,eventname):


  raw = bioread.read_file(filename[0])

  # print channel names
  for c, val in enumerate(list(raw.named_channels.keys())):
      print("Name of channel " + str(c + 1) + " is: " + val)

  channel_names = list(raw.named_channels.keys())
  for c, val in enumerate(channel_names):
      thefile = open(filename[:-5] + "_" + val + ".txt", 'w')
      res_dat = raw.channels[c].data[::int(raw.samples_per_second / dSF)]  # resample to 50 HZ
      for item in res_dat:
          thefile.write("%s\n" % item)
      thefile.close()

  # Get stimulus events
  eventchannel = [c for c, val in enumerate(channel_names) if val == eventname]
  stimcodes = raw.channels[eventchannel[0]].data
  stimcodes = [c for c, val in enumerate(np.diff(stimcodes)) if val == 5]
  stimcodes = [int(val / int(raw.samples_per_second / self.inputdata.Hz)) for val in stimcodes]
  thefile = open(filename[0][:-5] + "_stimcodes.txt", 'w')
  for item in stimcodes:
      thefile.write("%s\n" % item)
  thefile.close()

  print("files coverted :)")
