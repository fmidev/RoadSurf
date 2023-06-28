#pragma once

#include <boost/thread.hpp>

using MutexType = boost::shared_mutex;
using ReadLock = boost::shared_lock<MutexType>;
using WriteLock = boost::unique_lock<MutexType>;
using UpgradeReadLock = boost::upgrade_lock<MutexType>;
using UpgradeWriteLock = boost::upgrade_to_unique_lock<MutexType>;
