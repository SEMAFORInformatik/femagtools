# -*- coding: utf-8 -*-
"""
    femagtools.condor
    ~~~~~~~~~~~~~~~~~

    Creating and managing condor jobs



"""
#
import os
import subprocess
import sys
import datetime
import re
import time
import json
import fnmatch
import logging
import femagtools.job

logger = logging.getLogger(__name__)


def getDateTime(date, time):
    """complete date with year: m/d --> y/m/d"""
    if date == "???":
        return ""
    today = datetime.date.today()
    m, d = date.split('/')
    fulldate = today.replace(month=int(m), day=int(d))
    if fulldate > datetime.date.today():  # greater than today => decrease year
        fulldate = fulldate.replace(year=(datetime.date.today().year-1))
    return "{0}T{1}:00".format(fulldate.isoformat(), time)


class CondorCluster(object):
    """manages condor cluster directories"""

    def __init__(self, userdir):
        self.basedir = userdir
        self.dataFile = 'cluster.json'

    def getAvailableClusters(self):
        results = []
        for fn in self.getClusterDirectories():
            try:
                with open(fn, 'r') as f:
                    vars = json.load(f)
                results.append(vars['clusterId'])
            except ValueError:
                continue
        return results

    def getClusterDirectories(self):
        results = []

        for base, dirs, files in os.walk(self.basedir):
            goodfiles = fnmatch.filter(files, self.dataFile)
            for f in goodfiles:
                results.append(os.path.join(base, f))
        return results

    def getClusterData(self, directory):
        for d in self.getClusterDirectories():
            logger.info("Query %s != %s", os.path.basename(d), directory)
            if os.path.basename(os.path.dirname(d)) == directory:
                logger.info("Found %s", d)
                with open(d, 'r') as f:
                    vars = json.load(f)
                    return vars

        logger.info("Not Found %s", d)
        return None

    def getClusterList(self, clusterIds=[]):
        results = []
        for fn in self.getClusterDirectories():
            try:
                with open(fn, 'r') as f:
                    vars = json.load(f)
            except ValueError:
                continue

            # tooltip
            s = "<h3>ClusterId {0}</h3><table>\n".format(vars['clusterId'])
            attrs = dict(variable='Variable', baseValue='BaseValue',
                         bounds=['StartValue','EndValue'],
                         nSteps='nSteps',unit='Unit', desc='Desc')
            for r in vars['optimize']['decision_vars']:
                for a in attrs:
                    if isinstance(attrs[a], str):
                        rs = re.sub('@.*?@', '', str(r[a]))
                        s += "<tr><td>{0}:</td><td>{1}</td></tr>".format(attrs[a], rs)
                    elif isinstance(attrs[a], list):
                        for i in range(len(r[a])):
                            rs = re.sub('@.*?@', '', str(r[a][i]))
                            s += "<tr><td>{0}:</td><td>{1}</td></tr>".format(attrs[a][i], rs)
            s += "</table>"

            # decision_vars
            dv = []

            # get infos from each cluster directory
            if 'clusterId' in vars:
                obj = dict(clusterId=vars['clusterId'],
                           directory=os.path.basename(os.path.dirname(fn)),
                           timestamp=time.strftime("%Y-%m-%dT%H:%M:%S",
                                                   time.localtime(os.path.getctime(fn))),
                           tooltip=s,
                           decision_vars=vars['optimize']['decision_vars'],
                           selected=(vars['clusterId'] in clusterIds))
                results.append(obj)
        return results

    def status(self, clusterId=None):
        if clusterId:
            return self.getClusterData(clusterId)
        else:
            return self.getClusterList()


class Engine(object):
    """manages calculation tasks in a HTCondor environment"""

    def __init__(self):
        self.job = None

    def create_job(self, workdir):
        self.job = femagtools.job.CondorJob(workdir)
        return self.job

    def status(self, userdir=None):
        status_cmd = 'condor_status'
        try:
            proc = subprocess.Popen([status_cmd],
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)

            result = dict(machines=[], total=[])
            for l in proc.stdout:
                l = l.decode('utf-8').strip().split()
                if len(l) == 8 and \
                   (l[0].find(b'@') > 0 or str(l[5]).find('.') > 0):
                    result['machines'].append(dict(
                        name=l[0],
                        opsys=l[1],
                        arch=l[2],
                        state=l[3],
                        activity=l[4],
                        loadavg=float(l[5]),
                        mem=float(l[6]),
                        acttime=l[7]))
                elif len(l) == 8 and l[0] != b'Name':
                    result['total'].append(dict(oparch=l[0],
                                                total=int(l[1]),
                                                owner=int(l[2]),
                                                claimed=int(l[3]),
                                                unclaimed=int(l[4]),
                                                matched=int(l[5]),
                                                preempting=int(l[6]),
                                                backfill=int(l[7])))

            condorCl = CondorCluster(userdir)
            result['queue'] = self.queue(None)
            l = condorCl.getAvailableClusters()
            result['history'] = self.history(condorCl.getAvailableClusters())
            if userdir:
                result['cluster'] = condorCl.status()
            return result

        except ValueError as e:
            logger.error("Cannot invoke %s", status_cmd)
            e.args += status_cmd
            raise e

    def statusDir(self, userdir, directory):
        condorCl = CondorCluster(userdir)
        logger.info("%s, %s", userdir, directory)
        return condorCl.getClusterData(directory)

    def submit(self):
        submitfile = self.job.prepareDescription()
        cmdout = subprocess.check_output(["condor_submit", submitfile])
        self.clusterId = re.findall('\d+', cmdout.decode('utf-8'))[-1]
        logger.info('submit cluster %s', self.clusterId)
        return self.clusterId

    def join(self):
        """wait for all tasks to be terminated and return status"""
        ret = []
        if self.clusterId:
            cmd = ["condor_q", self.clusterId]
            while True:
                time.sleep(2)
                cmdout = subprocess.check_output(cmd)
                tasks = re.findall(self.clusterId+'\.\d+',
                                   cmdout.decode('utf-8'))
                if len(tasks) == 0:
                    break

            status = dict()
            cmdout = subprocess.check_output(
                ["condor_history", self.clusterId])
            for jobinfo in re.findall(
                    '^\s*{}\.\d+.+$'.format(
                        self.clusterId),
                    cmdout.decode('utf-8'), re.M):
                l = jobinfo.split()
                taskid = int(l[0].split('.')[-1])
                status[taskid] = l[5]
                logger.info('status %d: %s', taskid, l[5])
                self.job.setExitStatus(taskid, status[taskid])

            for k in sorted(status.keys()):
                ret.append(status[k])
            logger.info('finished cluster %s', self.clusterId)
        else:
            logger.warn('no condor cluster')

        return ret

    def queue(self, clusterId):
        cmd = ["condor_q"]
        results = []
        if clusterId:
            cmd.append(str(clusterId))
        proc = subprocess.Popen(cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)

        for l in proc.stdout:
            l = l.decode('utf-8').strip()
            if l.startswith('ID') or l.find('completed') > 0:
                continue
            props = l.split()
            if len(props) >= 9:
                s = props[4].split('+')
                t = s[1].split(':')
                c, p = props[0].split('.')
                results.append(dict(
                    id=props[0],
                    clusterId=c,
                    processId=p,
                    owner=props[1],
                    submitted=getDateTime(props[2], props[3]),
                    cpuUsage='{0:02}:{1}:{2}'.format(
                        int(24*s[0]+t[0]), t[1], t[2]),
                    status=props[5]))

        return results

    def history(self, clusterId):
        cmd = ["condor_history"]
        results = []
        if isinstance(clusterId, list):
            for clId in clusterId:
                results += self.history(str(clId))
            return results
        elif isinstance(clusterId, str):
            cmd.append(clusterId)

        logger.info("execute history cmd: " + str(cmd))
        proc = subprocess.Popen(cmd,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)

        for l in proc.stdout:
            l = l.decode('utf-8').strip()
            if l.startswith('ID') or l.find('completed') > 0:
                continue
            props = l.split()
            if len(props) >= 8:
                s = props[4].split('+')
                t = s[1].split(':')
                c, p = props[0].split('.')
                results.append(dict(
                    id=props[0],
                    clusterId=c,
                    processId=p,
                    owner=props[1],
                    submitted=getDateTime(props[2], props[3]),
                    cpuUsage='{0:02}:{1}:{2}'.format(
                        int(24*s[0] + t[0]), t[1], t[2]),
                    status=props[5],
                    completed=getDateTime(props[6], props[7])))
            else:
                sys.stderr.write("unknown condor_history line '{}'\n".format(
                    l))
        return results

