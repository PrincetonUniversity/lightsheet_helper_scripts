# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 15:00:21 2019

@author: lvbt
"""

import os, time, smtplib, ssl
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

#this is run on the ls acquisition computer

def send_email(img, receiver_email, check):
    """ sends completion email based on the time when the 50th to last plane was modified """
    #get stats
    f = os.stat(os.path.join(img, os.listdir(img)[len(os.listdir(img))-50]))
    
    if receiver_email == "": receiver_email = "zmd@princeton.edu"
    if check == "": check = 300
    else: check = int(check)    
    
    while time.time() - f.st_mtime < 180: #if its not been modified for over 3 minutes 
        
        print("\nlast modified time: ", time.ctime(f.st_mtime))
        
        print("\nacquisition not complete, will check again in {} minutes".format(round((check/60), 2)))
        time.sleep(check) #sleep for 5 minutes
        
        #check stats again
        f = os.stat(os.path.join(img, os.listdir(img)[len(os.listdir(img))-50]))
        
    #send email
    smtp_server = "smtp.gmail.com"
    port = 587  # For starttls
    sender_email = "pni.lvbt.lightsheet@gmail.com"
    password = "stephanisawesome"
    msg = MIMEMultipart()
    msg["Subject"] = "Light-sheet acquisition complete!"
    body = "Last modified time of file: {}".format(time.ctime(f.st_mtime))
    msg.attach(MIMEText(body, 'plain'))
    context = ssl.create_default_context()
    with smtplib.SMTP(smtp_server, port) as server:
        server.ehlo()  # Can be omitted
        server.starttls(context=context)
        server.ehlo()  # Can be omitted
        server.login(sender_email, password)
        text = msg.as_string()
        server.sendmail(sender_email, receiver_email, text)
    print ("\nsent email :D \n")
    

if __name__ == "__main__":
    
    img = "G:\\190308_20170212_tp_bl6_crii_1000r_02_4x_647_008na_1hfds_z7d5um_75msec_10povlp_08-31-27"
    
    print("press enter twice to accept all defaults when prompted for input")
    
    receiver_email = input("enter receiver email (default is zmd@princeton.edu): ")
    
    check = input("seconds after which you would like to check for updates (default is 300 s or 5 min): ")
    
    send_email(img, receiver_email, check)
    
    